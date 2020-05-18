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