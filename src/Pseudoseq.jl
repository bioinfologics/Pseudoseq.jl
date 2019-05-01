module Pseudoseq

export
    # High-level API
    
    ## Build-a-genome
    
    ### Chromosome blueprint
    plan_chrom,
    plan_repetition,
    plan_het,
    
    available_blocks,
    isavailable,
    consume_blocks,
    
    fabricate,
    
    ### Utility methods
    suggest_regions,
    suggest_alleles,
    count_het,
    
    # Sequencing
    
    ## High level user methods
    sequence,
    
    ## High level API
    ### Molecule Pool
    makepool,
    needed_sample_size,
    
    #### Transformations
    amplify,
    fragment,
    subsample,
    
    ### Reads
    make_reads,
    #### Technologies
    PairedEnd,
    SingleEnd,
    TaggedPairs,
    
    #### Information functions
    total_bp,
    
    #### Transformations
    mark_errors,
    mark_errors!,
    generate,
    
    ### SequencingView
    SequencingView,
    isfwd,
    isbwd,
    seqid,
    tag,
    iscircular,
    
    ### Views
    Views
    
using BioSequences
import Random

# Build-a-Genome
include("build_a_genome/internals.jl")
include("build_a_genome/chrom_blueprint.jl")
include("build_a_genome/plan/repetitions.jl")
include("build_a_genome/plan/polymorphism.jl")
include("build_a_genome/plan/chromosomes.jl")
include("build_a_genome/fabrication.jl")

# Operations on views
include("views/sequencing_view.jl")
include("views/views.jl")
include("views/transformations.jl")
include("views/summarize.jl")

needed_sample_size(C::Int, G::Int, L::Int) = div(C * G, L)
expected_coverage(G::Int, L::Int, N::Int) = div(L * N, G)

include("pool.jl")
include("reads.jl")
include("sequence.jl")

end # module
