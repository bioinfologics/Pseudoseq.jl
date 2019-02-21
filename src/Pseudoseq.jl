module Pseudoseq

export
    makeuniverse,
    fragment,
    sample,
    needed_sample_size,
    make_reads,
    make_paired_reads,
    make_errors,
    
    SequencingView,
    iscircular,
    isfwd,
    isbwd,
    seqid,
    tag,
    Views,
    save_views,
    load_views,
    multiply_views,
    tag_views,
    save_views_and_errs,
    load_views_and_errors,
    extract_sequence,
    add_errors!

using BioSequences
using JLD
using Distributions
using UnicodePlots
import Random.randperm

include("views/sequencing_view.jl")
include("views/views.jl")
include("views/multiply.jl")
include("views/sample.jl")
include("views/fragment.jl")

include("summarize.jl")

include("tag.jl")

include("errors.jl")
include("reads.jl")

include("universe.jl")




end # module
