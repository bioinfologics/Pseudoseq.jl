
struct BreakpointFragmenter
    seq::SequencingView
    nbreak::Int
end

Base.IteratorSize(::Type{BreakpointFragmenter}) = Base.HasLength()
Base.IteratorEltype(::Type{BreakpointFragmenter}) = Base.HasEltype()
Base.eltype(::Type{BreakpointFragmenter}) = SequencingView
Base.length(it::BreakpointFragmenter) = it.nbreak + 1

function Base.iterate(it::BreakpointFragmenter)
    breakpoints = sort!(resize!(randperm(length(it.seq)), it.nbreak))
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
    summarize_lengths(newviews)
    return newviews
end