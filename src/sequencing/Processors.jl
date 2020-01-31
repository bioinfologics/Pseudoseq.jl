


struct Amplifier
    n::Int
end
amplify(n::Int) = Amplifier(n)
(a::Amplifier)(p::Molecules) = amplify(p, a.n)

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
(s::SubSampler)(p::Molecules) = subsample(p, s.n)

struct Selector{F<:Function}
    f::F
end
select(f::Function) = Selector(f)
(s::Selector{F})(p::Molecules) where {F<:Function} = select(s.f, p)