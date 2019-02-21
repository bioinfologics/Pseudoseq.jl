
struct SequencingView
    id::UInt64
    start::UInt32
    stop::UInt32
    tag::UInt64
end

function SequencingView(id::Integer, start::Integer, stop::Integer, circular::Bool = false)
    newid = (convert(UInt64, id) << 1) | convert(UInt64, circular)
    return SequencingView(newid, convert(UInt32, start), convert(UInt32, stop), zero(UInt64))
end

@inline iscircular(x::SequencingView) = convert(Bool, x.id & one(typeof(x.id)))
@inline seqid(x::SequencingView) = x.id >> 1
@inline tag(x::SequencingView) = x.tag
@inline window(x::SequencingView) = x.start:x.stop
@inline isfwd(x::SequencingView) = x.start < x.stop
@inline isbwd(x::SequencingView) = x.start > x.stop

@inline Base.length(x::SequencingView) = ifelse(isfwd(x), x.stop - x.start + 1, x.start - x.stop + 1)
@inline Base.eachindex(x::SequencingView) = Base.OneTo(length(x))
@inline Base.firstindex(x::SequencingView) = 1
@inline Base.lastindex(x::SequencingView) = length(x)
@inline Base.first(x::SequencingView) = x.start
@inline Base.last(x::SequencingView) = x.stop

function Base.show(io::IO, x::SequencingView)
    print(io, "A ", length(x), "bp view over a sequence ")
    if isbwd(x)
        print(io, '(', last(x), " <=== ", first(x), ')')
    else
        print(io, '(', first(x), " ===> ", last(x), ')')
    end
end
    
@inline function Base.checkbounds(x::SequencingView, i::Integer)
    if 1 ≤ i ≤ length(x)
        return true
    end
    throw(BoundsError(x, i))
end

@inline function Base.getindex(x::SequencingView, i::Integer)
    @boundscheck checkbounds(x, i)
    return ifelse(isfwd(x), first(x) + i - 1, first(x) - i + 1)
end

function subview(x::SequencingView, start::Integer, stop::Integer)
    return SequencingView(x.id, x[start], x[stop], tag(x))
end

Base.reverse(x::SequencingView) = SequencingView(seqid(x), last(x), first(x), tag(x))

    
    
    