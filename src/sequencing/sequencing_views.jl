#
# Sequencing views
# ================


# Abstract type and high level behaviour
# --------------------------------------

abstract type AbstractSequencingView end

# Specific concrete view types
# ----------------------------

# A basic view type - used most often.

struct BasicSequencingView <: AbstractSequencingView
    sequence::UInt64
    start::UInt32
    stop::UInt32
end

function BasicSequencingView(id::Integer, start::Integer, stop::Integer)
    return BasicSequencingView(convert(UInt64, id), convert(UInt32, start), convert(UInt32, stop))
end

@inline function BasicSequencingView(bsv::BasicSequencingView, start::Integer, stop::Integer)
    return BasicSequencingView(seqid(bsv), start, stop)
end

@inline seqid(x::BasicSequencingView) = x.sequence
@inline Base.first(x::BasicSequencingView) = x.start
@inline Base.last(x::BasicSequencingView) = x.stop

# A view type that includes meta info in the form of an additional UInt64 tag.

struct TaggedSequencingView <: AbstractSequencingView
    sequence::UInt64
    start::UInt32
    stop::UInt32
    tag::UInt64
end

function TaggedSequencingView(id::Integer, start::Integer, stop::Integer, tag::UInt64 = zero(UInt64))
    return TaggedSequencingView(id, convert(UInt32, start), convert(UInt32, stop), tag)
end

@inline function TaggedSequencingView(tsv::TaggedSequencingView, start::Integer, stop::Integer)
    return TaggedSequencingView(seqid(tsv), start, stop, tag(tsv))
end

@inline seqid(x::TaggedSequencingView) = x.sequence
@inline Base.first(x::TaggedSequencingView) = x.start
@inline Base.last(x::TaggedSequencingView) = x.stop
@inline tag(x::TaggedSequencingView) = x.tag

const BasicViews = Vector{BasicSequencingView}
const TaggedViews = Vector{TaggedSequencingView}

# Generic sequencing view methods
# -------------------------------

@inline isfwd(x::AbstractSequencingView) = first(x) < last(x)
@inline isbwd(x::AbstractSequencingView) = first(x) > last(x)

@inline function Base.length(x::AbstractSequencingView)
    return ifelse(isfwd(x), last(x) - first(x) + 1, first(x) - last(x) + 1)
end

@inline Base.firstindex(x::AbstractSequencingView) = 1
@inline Base.lastindex(x::AbstractSequencingView) = length(x)
@inline window(x::AbstractSequencingView) = first(x):last(x)
@inline Base.eachindex(x::AbstractSequencingView) = Base.OneTo(lastindex(x))

function Base.show(io::IO, x::AbstractSequencingView)
    print(io, "A ", length(x), "bp view over a sequence ")
    if isbwd(x)
        print(io, '(', last(x), " <=== ", first(x), ')')
    else
        print(io, '(', first(x), " ===> ", last(x), ')')
    end
end

@inline function Base.checkbounds(x::AbstractSequencingView, i::Integer)
    if 1 ≤ i ≤ length(x)
        return true
    end
    throw(BoundsError(x, i))
end

@inline function Base.getindex(x::AbstractSequencingView, i::Integer)
    @boundscheck checkbounds(x, i)
    return ifelse(isfwd(x), first(x) + i - 1, first(x) - i + 1)
end

@inline function Base.reverse(x::T) where {T<:AbstractSequencingView}
    return T(x, last(x), first(x))
end

@inline function subview(x::T, start::Integer, stop::Integer) where {T<:AbstractSequencingView}
    return T(x, x[start], x[stop])
end

function take_paired_ends(sv::AbstractSequencingView, flen::Int, rlen::Int)
    return subview(sv, firstindex(sv), flen), subview(sv, lastindex(sv), lastindex(sv) - rlen + 1)
end

function take_single_end(sv::AbstractSequencingView, len::Int = length(sv))
    if (len > length(sv)) | (len < 0)
        throw(ArgumentError("len must be 0 < len ≤ length(sv)"))
    end
    forwards = rand(Bool)
    start = ifelse(forwards, firstindex(sv), lastindex(sv))
    stop = ifelse(forwards, len, lastindex(sv) - len + 1)
    return subview(sv, start, stop)
end

function extract_sequence(ref::BioSequence, view::AbstractSequencingView)
    if isbwd(view)
        subseq = ref[last(view):first(view)]
        reverse_complement!(subseq)
    else
        subseq = ref[first(view):last(view)]
    end
    return subseq
end

function extract_sequence(genome::Vector{<:BioSequence}, view::AbstractSequencingView)
    return extract_sequence(genome[seqid(view)], view)
end

# Utils for dealing with multiple views
# -------------------------------------

# Transformations

## The amplify transformation

function amplify(views::Vector{T}, n::Integer) where {T<:AbstractSequencingView}
    newviews = Vector{T}(undef, length(views) * n)
    i = 0
    for view in views
        @inbounds for _ in 1:n
            i += 1
            newviews[i] = view
        end
    end
    return newviews
end

# Utility function, randomly select values from a range, without replacement.
function sample_values(range::UnitRange, n::Int)
    @assert n <= length(range)
    s = Set{eltype(range)}()
    while length(s) < n
        push!(s, rand(range))
    end
    return collect(s)
end

## The subsample transformation

function subsample(vs::Vector{T}, n::Int) where {T<:AbstractSequencingView}
    if length(vs) < n
        throw(ArgumentError("number of views to sample is greater than total number of views."))
    end
    p = sample_values(1:length(vs), n)
    return vs[p]
end

## The selection transformation

function select(fn::Function, vs::Vector{T}) where {T<:AbstractSequencingView}
    nselected = 0
    newviews = Vector{T}(undef, length(vs))
    @inbounds for oldview in vs
        prob = fn(oldview)
        if rand() < prob
            nselected += 1
            newviews[nselected] = oldview
        end
    end
    resize!(newviews, nselected)
    return newviews
end

## The tag transformation

function tag(vs::BasicViews, npools::Int)
    newviews = TaggedViews(undef, length(vs))
    tags = rand(one(UInt64):UInt64(npools), length(vs))
    for i in eachindex(vs)
        view = vs[i]
        newviews[i] = TaggedSequencingView(seqid(view), first(view), last(view), tags[i])
    end
    return newviews
end

## The fragment transformation

struct BreakpointFragmenter{T<:AbstractSequencingView}
    seq::T
    nbreak::Int
end

Base.IteratorSize(::Type{BreakpointFragmenter{T}}) where {T} = Base.HasLength()
Base.IteratorEltype(::Type{BreakpointFragmenter{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{BreakpointFragmenter{T}}) where T = T
Base.length(it::BreakpointFragmenter) = it.nbreak + 1

function Base.iterate(it::BreakpointFragmenter{T}) where {T}
    breakpoints = sort!(sample_values(1:length(it.seq), it.nbreak))
    return iterate(it, (1, 0, breakpoints))
end

function Base.iterate(it::BreakpointFragmenter{T}, state::Tuple{Int,Int,Vector{Int}}) where {T}
    i, lastbp, breakpoints = state
    lastbp == lastindex(it.seq) && return nothing
    newstart = lastbp + 1
    newstop = i > it.nbreak ? lastindex(it.seq) : breakpoints[i]
    return subview(it.seq, newstart, newstop), (i + 1, newstop, breakpoints)
end

fragments(seq::AbstractSequencingView, nbreak::Int) = BreakpointFragmenter{typeof(seq)}(seq, nbreak)

function fragment(vs::Vector{T}, meansize::Int) where {T<:AbstractSequencingView}
    newviews = Vector{T}()
    for view in vs
        nbreak = div(length(view), meansize) - 1
        for f in fragments(view, nbreak)
            push!(newviews, f)
        end
    end
    return newviews
end

## End-taking transformations

function take_paired_ends(vs::Vector{T}, flen::Int, rlen::Int) where {T<:AbstractSequencingView}
    i = 0
    newviews = Vector{T}(undef, length(vs) * 2)
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

function take_single_ends(vs::Vector{T}, len::Int) where {T<:AbstractSequencingView}
    newviews = Vector{T}(undef, length(vs))
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

function take_single_ends(vs::Vector{T}, len::Nothing) where {T<:AbstractSequencingView}
    newviews = Vector{T}(undef, length(vs))
    for (i, view) in enumerate(vs)
        newviews[i] = take_single_end(view)
    end
    return newviews
end

## Flipping transformation
#=
function flip_views!(vs::Vector{<:AbstractSequencingView})
    willflip = rand(Bool, length(vs))
    @inbounds for (i, view) in enumerate(vs)
        if willflip[i]
            vs[i] = reverse(vs[i])
        end
    end
    return vs
end

flip_views(vs::Vector{<:AbstractSequencingView}) = flip_views!(deepcopy(vs))
=#
# Summarizing functions
# ---------------------

function summarize_lengths(vs::Vector{<:AbstractSequencingView})::Tuple{Int, Int, Int}
    lens = length.(vs)
    return maximum(lens), div(sum(lens), length(vs)), minimum(lens)
end

function summarize_tags(vs::Vector{TaggedSequencingView})
    tagdict = Dict{UInt64, Int}()
    for v in vs
        t = tag(v)
        if t > 0
            tagdict[t] = get(tagdict, t, 0) + 1
        end
    end
    return tagdict
end
