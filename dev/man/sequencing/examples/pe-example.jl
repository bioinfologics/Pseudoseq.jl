using Pseudoseq.Sequencing

sequence("ecoli-ref.fasta", "pe-reads.fastq"; ng = 5000, flen = 700, cov = 50, paired = true, rdlen = 250, err = 0.001)

pool = Molecules("ecoli-ref.fasta", 5000)

cutpool = fragment(pool, 700)

genome_size = 4639675
expected_coverage = 50
read_length = 250

N = needed_sample_size(expected_coverage, genome_size, read_length)
N = div(N, 2) # Divide by 2 as we're doing paired end sequencing.

sampledpool = subsample(cutpool, N)

pe_reads = paired_reads(sampledpool, 250)

pe_w_errs = mark_errors(pe_reads, 0.001)

generate("pe-reads.fastq", pe_w_errs)

pool = Molecules("ecoli-ref.fasta", 5000)

cutter = fragment(700)
sampler = subsample(N) # Remember how to computed N previously.
mkreads = paired_reads(250)
adderr = mark_errors(0.001)

pool |> cutter |> sampler |> mkreads |> adderr |> generate("pe-reads.fastq")

my_protocol = adderr ∘ mkreads ∘ sampler ∘ cutter

pool |> my_protocol |> generate("pe-reads.fastq")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

