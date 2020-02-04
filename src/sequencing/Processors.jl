# Processors
# ==========

# Processors for pipelining


struct Amplifier{N}
    n::NTuple{N,Int}
end
amplify(n::Int...) = Amplifier{length(n)}(n)
(a::Amplifier)(p::Molecules) = amplify(p, a.n)
function (a::Amplifier{N})(p::Molecules...) where {N}
    return tuple(amplify(sub_p, n_sub) for n_sub in a.n for sub_p in p)...)
end

struct Fragmenter
    meansize::Int
end
fragment(meansize::Int) = Fragmenter(meansize)
(f::Fragmenter)(p::Molecules) = fragment(p, f.meansize)

struct Tagger
    ntags::Int
end
tag(ntags::Int) = Tagger(ntags)
(t::Tagger)(p::Molecules) = tag(p, t.ntags)

struct SubSampler
    nsamples::Int
end
subsample(n::Int) = SubSampler(n)
(s::SubSampler)(p::Molecules) = subsample(p, s.nsamples)

struct Selector{F<:Function}
    f::F
end
select(f::Function) = Selector(f)
(s::Selector{F})(p::Molecules) where {F<:Function} = select(s.f, p)

struct PairedReads
    flen::Int
    rlen::Int
end
paired_reads(flen::Int, rlen::Int = flen) = PairedReads(flen, rlen)
(pr::PairedReads)(p::Molecules) = paired_reads(p, pr.flen, pr.rlen)

struct UnPairedReads{T<:Union{Nothing,Int}}
    len::T
end
unpaired_reads(len::T) where {T<:Union{Nothing,Int}} = UnPairedReads(len)
(ur::UnPairedReads)(p::Molecules) = unpaired_reads(p, ur.len)

struct ErrorMaker
    rate::Float64
end
mark_errors(rate::Float64) = ErrorMaker(rate)
(em::ErrorMaker)(p::Reads) = mark_errors(p, em.rate)

struct FileGenerator
    filename::String
end
generate(filename::String) = FileGenerator(filename)
(fg::FileGenerator)(r::Reads) = generate(fg.filename, r)