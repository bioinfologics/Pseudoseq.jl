# The amplify transformation
# ---------------------------

function amplify(views::Views, n::Integer)
    newviews = Views(length(views) * n)
    i = 0
    for view in views
        @inbounds for _ in 1:n
            i += 1
            newviews[i] = view
        end
    end
    return newviews
end


# The subsample transformation
# ----------------------------

function subsample(vs::Views, n::Int)
    if length(vs) < n
        throw(ArgumentError("number of views to sample is greater than total number of views."))
    end
    p = sample_values(1:length(vs), n)
    return vs[p]
end


# The tag transformation
# ----------------------

function tag(vs::Views, npools::Int)
    newviews = Views(length(vs))
    tags = rand(one(UInt64):UInt64(npools), length(vs))
    for i in eachindex(vs)
        view = vs[i]
        newviews[i] = SequencingView(view.id, first(view), last(view), tags[i])
    end
    return newviews
end


# The fragment transformation
# ---------------------------

struct BreakpointFragmenter
    seq::SequencingView
    nbreak::Int
end

Base.IteratorSize(::Type{BreakpointFragmenter}) = Base.HasLength()
Base.IteratorEltype(::Type{BreakpointFragmenter}) = Base.HasEltype()
Base.eltype(::Type{BreakpointFragmenter}) = SequencingView
Base.length(it::BreakpointFragmenter) = it.nbreak + 1

function sample_values(range::UnitRange, n::Int)
    @assert n <= length(range)
    s = Set{eltype(range)}()
    while length(s) < n
        push!(s, rand(range))
    end
    return collect(s)
end

function Base.iterate(it::BreakpointFragmenter)
    breakpoints = sort!(sample_values(1:length(it.seq), it.nbreak))
    return iterate(it, (1, 0, breakpoints))
end

function Base.iterate(it::BreakpointFragmenter, state::Tuple{Int,Int,Vector{Int}})
    i, lastbp, breakpoints = state
    lastbp == lastindex(it.seq) && return nothing
    newstart = lastbp + 1
    newstop = i > it.nbreak ? lastindex(it.seq) : breakpoints[i]
    return subview(it.seq, newstart, newstop), (i + 1, newstop, breakpoints)
end

fragments(seq::SequencingView, nbreak::Int) = BreakpointFragmenter(seq, nbreak)

function fragment(vs::Views, meansize::Int)
    newviews = Views()
    for view in vs
        nbreak = div(length(view), meansize) - 1
        for f in fragments(view, nbreak)
            push!(newviews, f)
        end
    end
    return newviews
end

# End-taking transformations
# --------------------------

# These are basically just used to prepare a Views struct for creating a Reads struct.

function take_paired_ends(vs::Views, flen::Int, rlen::Int)
    i = 0
    newviews = Views(length(vs) * 2)
    for view in vs
        vl = length(view)
        if vl > flen && vl > rlen # Guard against the sloppy fragmenter.
            fwdview, revview = take_paired_ends(view, flen, rlen)
            newviews[i += 1] = fwdview
            newviews[i += 1] = revview
        end
    end
    resize!(newviews, i)
    return newviews
end

function take_single_ends(vs::Views, len::Int)
    newviews = Views(length(vs))
    i = 0
    for view in vs
        vl = length(view)
        if vl > len
            newviews[i += 1] = take_single_end(view, len)
        end
    end
    resize!(newviews, i)
    return newviews
end

function take_single_ends(vs::Views)
    newviews = Views(length(vs))
    for (i, view) in enumerate(vs)
        newviews[i] = Pseudoseq.take_single_end(view)
    end
    return newviews
end