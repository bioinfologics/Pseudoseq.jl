using Pseudoseq

sequence("ecoli-ref.fasta", "longreads.fastq"; ng = 5000, flen = 40000, cov = 30, rdlen = nothing, err = 0.1, paired = false)

pool = makepool("ecoli-ref.fasta", 5000)

cutpool = fragment(pool, 40000)

genome_size = 4639675
expected_coverage = 30
readlength = 40000

N = needed_sample_size(expected_coverage, genome_size, readlength)

sampledpool = subsample(cutpool, N)

se_reads = make_reads(SingleEnd, sampledpool)

se_w_errs = mark_errors(se_reads, 0.1)

generate("longreads.fastq", se_w_errs)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

