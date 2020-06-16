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
    makereads


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

function needed_sample_size(coverage::ExpectedCoverage, genome_len::SequenceLength, lens)
    return div(coverage * genome_len, sum(lens))
end



expected_coverage(G::Int, L::Int, N::Int) = div(L * N, G)

include("sequencing_views.jl")
include("Molecules.jl")
include("Reads.jl")
include("DSL.jl")
include("sequence.jl")

end # module