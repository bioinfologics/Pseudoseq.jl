
struct Views
    viewid::Vector{UInt64}
    viewstart::Vector{UInt32}
    viewstop::Vector{UInt32}
    tags::Vector{UInt64}
end

Views() = Views(Vector{UInt64}(), Vector{UInt32}(), Vector{UInt32}(), Vector{UInt64}())

function Views(n::Integer)
    return Views(Vector{UInt64}(undef, n), Vector{UInt32}(undef, n), Vector{UInt32}(undef, n), Vector{UInt64}(undef, n))
end

# TODO: May not be stable across versions of BioSequences.jl or FASTX.jl
function Views(input::FASTA.Reader, circular::Bool = false)
    views = Views()
    i = 1
    rec = eltype(input)()
    while !eof(input)
        read!(input, rec)
        push!(views, SequencingView(i, 1, length(rec.sequence), circular))
        i += 1
    end
    return views
end

@inline function Base.getindex(vs::Views, i::Integer)
    return SequencingView(vs.viewid[i], vs.viewstart[i], vs.viewstop[i], vs.tags[i])
end

function Base.getindex(vs::Views, is::Vector{<:Integer})
    newviews = Views(length(is))
    for (i, j) in enumerate(is)
        newviews[i] = vs[j]
    end
    return newviews
end

@inline function Base.setindex!(vs::Views, sv::SequencingView, i::Integer)
    vs.viewid[i] = sv.id
    vs.viewstart[i] = first(sv)
    vs.viewstop[i] = last(sv)
    vs.tags[i] = tag(sv)
    return sv
end

@inline function Base.iterate(vs::Views, state = 1)
    if state > length(vs)
        nothing
    else
        newstate = state + 1
        return vs[state], newstate
    end
end

@inline Base.length(vs::Views) = length(vs.viewid)
@inline Base.eachindex(vs::Views) = Base.OneTo(length(vs))

function Base.push!(x::Views, y::SequencingView)
    push!(x.viewid, y.id)
    push!(x.viewstart, first(y))
    push!(x.viewstop, last(y))
    push!(x.tags, tag(y))
    return y
end

function Base.resize!(x::Views, s::Integer)
    resize!(x.viewid, s)
    resize!(x.viewstart, s)
    resize!(x.viewstop, s)
    resize!(x.tags, s)
end

function save_views(vs::Views, prefix::String)
    ds = jldopen(string(first(splitext(prefix)), ".jld"), "w")
    write_views(ds, vs)
    close(ds)
end

function save_views_and_errs(vs::Views, errs::Tuple{Vector{Int},Vector{Int}}, prefix::String)
    ds = jldopen(string(first(splitext(prefix)), ".jld"), "w")
    write_views(ds, vs)
    write(ds, "errread", errs[1])
    write(ds, "errpos", errs[2])
    close(ds)
end

function write_views(ds::JLD.JldFile, vs::Views)
    write(ds, "ids", vs.viewid)
    write(ds, "starts", vs.viewstart)
    write(ds, "stops", vs.viewstop)
    write(ds, "tags", vs.tags)
end

function load_views(file::String)
    ds = jldopen(file, "r")
    ids = read(ds, "ids")
    starts = read(ds, "starts")
    stops = read(ds, "stops")
    tags = read(ds, "tags")
    close(ds)
    return Views(ids, starts, stops, tags)
end

function load_views_and_errors(file::String)
    ds = jldopen(file, "r")
    vs = read_views(ds)
    errread = read(ds, "errread")
    errpos = read(ds, "errpos")
    close(ds)
    
    # Populate error dictionary
    errs = Dict{Int, Vector{Int}}()
    for (read, pos) in zip(errread, errpos)
        if haskey(errs, read)
            push!(errs[read], pos)
        else
            errs[read] = [pos]
        end
    end
    return vs, errs
end

function read_views(ds::JLD.JldFile)
    ids = read(ds, "ids")
    starts = read(ds, "starts")
    stops = read(ds, "stops")
    tags = read(ds, "tags")
    return Views(ids, starts, stops, tags)
end
