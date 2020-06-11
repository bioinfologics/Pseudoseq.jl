# Processors
# ==========

# Processors for pipelining


struct Amplifier
    n::Int
end
amplify(n::Int) = Amplifier(n)
(a::Amplifier)(p::Molecules) = amplify(p, a.n)

struct Fragmenter
    meansize::Int
end
fragment(meansize::SequenceLength) = Fragmenter(meansize.val)
(f::Fragmenter)(p::Molecules) = fragment(p, f.meansize)

struct Tagger
    ntags::Int
end
tag(ntags::Int) = Tagger(ntags)
(t::Tagger)(p::Molecules) = tag(p, t.ntags)

struct NSubSampler
    nsamples::Int
end
subsample(n::Int) = NSubSampler(n)
(s::NSubSampler)(p::Molecules) = subsample(p, s.nsamples)

struct CovSubSampler{T<:Union{SequenceLength,Tuple{SequenceLength,SequenceLength}}}
    cov::ExpectedCoverage
    readlen::T
end
subsample(cov::ExpectedCoverage, rlens)  = CovSubSampler{typeof(rlens)}(cov, rlens)
(s::CovSubSampler{T})(p::Molecules) where {T} = subsample(p, s.cov, s.readlen)

struct Selector{F<:Function}
    f::F
end
select(f::Function) = Selector(f)
(s::Selector{F})(p::Molecules) where {F<:Function} = select(s.f, p)

makereads(len::SequenceLength) = unpaired_reads(len.val)
makereads() = unpaired_reads(nothing)
makereads(lena::SequenceLength, lenb::SequenceLength) = paired_reads(lena.val, lenb.val)
makereads(lens::Tuple{SequenceLength,SequenceLength}) = makereads(lens...)


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

struct SubstitutionMaker{F<:Function}
    fun::F
end
make_substitutions(f::F) where {F<:Function} = SubstitutionMaker{F}(f)
(sm::SubstitutionMaker{F})(p::Reads) where {F<:Function} = edit_substitutions(sm.fun, p)

struct FileGenerator
    filename::String
end
generate(filename::String) = FileGenerator(filename)
(fg::FileGenerator)(r::Reads) = generate(fg.filename, r)