using Pseudoseq

dnapool = makepool("../ecoli-ref.fasta", 5000)
cutpool = fragment(dnapool, 40000)

taggedpool = tag(cutpool, 1000000)

taggedcutpool = fragment(taggedpool, 700)

genome_size = 4639675
expected_coverage = 50
read_length = 250

N = needed_sample_size(expected_coverage, genome_size, read_length)
N = div(N, 2) # Divide by 2 as we're doing paired end sequencing.

sampledpool = subsample(taggedcutpool, N)

tagged_reads = make_reads(TaggedPairs, sampledpool, 250)
tagged_w_errs = mark_errors(tagged_reads, 0.001)

generate("tagged_reads.fastq", tagged_w_errs)#-
# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

