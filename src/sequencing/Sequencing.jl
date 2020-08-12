module Sequencing

export
    # Molecules
    Molecules,
    amplify,
    fragment,
    select,
    subsample,
    tag,
    flip,
    # Reads
    paired_reads,
    unpaired_reads,
    nreads,
    ClearSubstitutions,
    FixedProbSubstitutions,
    edit_substitutions,
    make_substitutions,
    generate,
    # Utils
    needed_sample_size,
    expected_coverage,
    sequence,
    X, x, bp,
    makereads,
    
    coverage_report


using BioSequences, FASTX, Random

struct ExpectedCoverage
    val::Int
end

const X = ExpectedCoverage
const x = ExpectedCoverage

Base.:*(val::Int, ::Type{ExpectedCoverage}) = ExpectedCoverage(val)
Base.show(io::IO, v::ExpectedCoverage) = println(io, v.val, "X sequencing coverage")

struct SequenceLength
    val::Int
end

const bp = SequenceLength

Base.:*(val::Int, ::Type{SequenceLength}) = SequenceLength(val)
Base.show(io::IO, x::SequenceLength) = println(io, x.val, " base pairs")

Base.:+(x::SequenceLength, y::SequenceLength) = SequenceLength(x.val + y.val)

Base.:*(cov::ExpectedCoverage, len::SequenceLength) = SequenceLength(cov.val * len.val)

function Base.:*(val::Int, rl::SequenceLength)
    @assert val == 2
    return (rl, rl)
end

Base.div(x::SequenceLength, y::SequenceLength) = div(x.val, y.val)
Base.iterate(x::SequenceLength) = (x, nothing)
Base.iterate(x::SequenceLength, ::Any) = nothing

function needed_sample_size(coverage::ExpectedCoverage, genome_len::SequenceLength, lens)
    return div(coverage * genome_len, sum(lens))
end



expected_coverage(G::Int, L::Int, N::Int) = div(L * N, G)

include("sequencing_views.jl")
include("Molecules.jl")
include("Reads.jl")
include("DSL.jl")
#include("sequence.jl")

struct CoverageReport
    genome::Vector{LongSequence{DNAAlphabet{2}}}
    covs::Vector{Vector{UInt}}
end

function coverage_report(x)
    covs = [zeros(UInt, length(g)) for g in genome(x)]
    vs = views(x)
    for v in vs
        c = covs[seqid(v)]
        @inbounds for i in min(first(v), last(v)):max(first(v), last(v))
            c[i] = c[i] + 1
        end
    end
    return CoverageReport(genome(x), covs)
end

function uncovered_positions(cr::Sequencing.CoverageReport, thresh = 0)
    return [[i for i in eachindex(c) if c[i] <= thresh] for c in cr.covs]
end

function _summarize_cov(V::Vector{UInt64})
    min = typemax(UInt64)
    max = typemin(UInt64)
    sum = typemin(UInt64)
    for v in V
        min = ifelse(v < min, v, min)
        max = ifelse(v > max, v, max)
        sum = sum + v
    end
    return (min, max, sum)
end

function summarize(cr::Sequencing.CoverageReport, bychrom = false)
    stats = [_summarize_cov(c) for c in cr.covs]
    if bychrom
        println("Coverage Summary:")
        for i in eachindex(stats)
            min, max, total = stats[i]
            mean = total / length(cr.covs[i])
            println("Chromosome $i:")
            println("\tmin: $min")
            println("\tmean: $mean")
            println("\tmax: $max")
        end
    else
        min = minimum(x[1] for x in stats)
        total = sum(x[3] for x in stats)
        max = maximum(x[2] for x in stats)
        mean = total / sum(length(x) for x in cr.covs)
        println("Coverage Summary:")
        println("\tmin: $min")
        println("\tmean: $mean")
        println("\tmax: $max")
    end
end

end # module