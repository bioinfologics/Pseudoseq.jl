"""
    sequence(input, output = nothing; kwargs...)
    
Run a sequencing experiment from start to finish.

This method is the highest-level function for running sequencing simulations
with Pseudoseq. It runs the Pseudoseq sequencin worksflow described in the manual,
tuned by the following keyword parameters:

Global parameters:

- `input`: A filename or `FASTA.Reader` providing the input genome.
- `output`: A filename or `FASTA.Writer` providing a destination for the output reads.

Parameters for workflow step 1 (Create a pool of DNA molecules):
- `ng`: Integer; the number of genomes the molecule pool is initialized with.

Parameters for step 2 (Processing the DNA molecule pool):
- `tusize`: Integer; the number of possible tags a DNA molecule may be tagged with (default: 0).
- `taggedflen`: Integer; the desired average fragment length at which molecules are tagged (default: 0bp). 
- `flen`: Integer; the desired average fragment length at which molecules are sampled (default: 700bp).
- `cov`: Integer; the desired expected coverage (default: 30x coverage).

Parameters for step 3 (Generating reads):
- `paired`: true or false; Whether or not to sequence from both end of each molecule (default: true).
- `rdlen`: Integer; the desired read length (default: 250bp).
- `err`: Float; the desired per-base error rate (default: 0.001).
"""
function sequence(input, output = nothing;
                  ng::Int = 1000,
                  tusize::Int = 0,
                  taggedflen::Int = 0,
                  flen::Int = 700,
                  cov::Int = 30,
                  paired::Bool = true,
                  rdlen = 250,
                  err::Float64 = 0.001
                  )
                  
    pool = makepool(input, 1)
    genome_size = sum(length(x) for x in views(pool))
    
    pool = amplify(pool, ng)
    
    @info string("- ✔ Created pool of ", ng, " copies of a ", genome_size, "bp genome")
    
    if taggedflen > 0 && tusize > 0
        pool = fragment(pool, taggedflen)
        @info string("- ✔ Created pool of fragments with an average length of ", taggedflen, "bp")
        pool = tag(pool, tusize)
        @info string("- ✔ Created pool of tagged fragments with an average length of ", taggedflen, "bp")
    end
    
    pool = fragment(pool, flen)
    @info string("- ✔ Created pool of fragments with an average length of ", flen, "bp")
    
    N = needed_sample_size(cov, genome_size, ifelse(rdlen === nothing, flen, rdlen))
    # Divide by 2 if we're doing paired end sequencing.
    N = ifelse(paired, div(N, 2), N)
    
    pool = subsample(pool, N)
    @info string("- ✔ Subsampled pool at ", cov, "X coverage (", N, " molecules)")
    
    if paired
        reads = make_reads(PairedEnd, pool, rdlen)
        @info string("- ✔ Created set of ", rdlen, "bp paired-end reads")
    else
        reads = make_reads(SingleEnd, pool, rdlen)
        @info string("- ✔ Created set of ", rdlen === nothing ? "" : string(rdlen, "bp "), "single-end reads")
    end
    
    reads_w_errs = mark_errors(reads, err)
    @info string("- ✔ Applied sequencing errors at a per-base rate of ", err)
    
    if output !== nothing
        generate(output, reads_w_errs)
    end
    
    return reads_w_errs
end