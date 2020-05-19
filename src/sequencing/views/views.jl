
struct Views
    viewstore::Vector{SequencingView}
end

Views() = Views(Vector{SequencingView}())

function Views(n::Integer)
    return Views(Vector{SequencingView}(undef, n))
end

@inline Base.getindex(vs::Views, i::Integer) = getindex(vs.viewstore, i)

function Base.getindex(vs::Views, is::Vector{<:Integer})
    newviews = Views(length(is))
    for (i, j) in enumerate(is)
        newviews[i] = vs[j]
    end
    return newviews
end

@inline Base.setindex!(vs::Views, sv::SequencingView, i::Integer) = setindex!(vs.viewstore, sv, i)

@inline function Base.iterate(vs::Views, state = 1)
    if state > length(vs)
        nothing
    else
        newstate = state + 1
        return vs[state], newstate
    end
end

@inline Base.length(vs::Views) = length(vs.viewstore)
@inline Base.eachindex(vs::Views) = Base.OneTo(length(vs))



Base.push!(x::Views, y::SequencingView) = push!(x.viewstore, y)

Base.resize!(x::Views, s::Integer) = resize!(x.viewstore, s)

Base.show(io::IO, x::Views) = show(io, "text/plain", x.viewstore)