module Sequencing

export
    Molecules,
    amplify,
    fragment,
    select,
    subsample,
    tag,
    flip,
    paired_reads,
    unpaired_reads,
    mark_errors,
    generate,
    needed_sample_size,
    expected_coverage,
    nreads,
    sequence
    

using BioSequences, FASTX, Random

needed_sample_size(C::Int, G::Int, L::Int) = div(C * G, L)
expected_coverage(G::Int, L::Int, N::Int) = div(L * N, G)

include("sequencing_views.jl")
include("Molecules.jl")
include("Reads.jl")
include("Processors.jl")
include("sequence.jl")

end # module