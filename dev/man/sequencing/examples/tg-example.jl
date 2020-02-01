using Pseudoseq.Sequencing

sequence("ecoli-ref.fasta", "tagged_reads.fastq"; ng = 5000, tusize = 1000000, taggedflen = 40000, flen = 700, cov = 50, rdlen = 250, err = 0.1)

dnapool = Molecules("ecoli-ref.fasta", 5000)

cutpool = fragment(dnapool, 40000)

taggedpool = tag(cutpool, 1000000)

taggedcutpool = fragment(taggedpool, 700)

genome_size = 4639675
expected_coverage = 50
read_length = 250

N = needed_sample_size(expected_coverage, genome_size, read_length)
N = div(N, 2) # Divide by 2 as we're doing paired end sequencing.

sampledpool = subsample(taggedcutpool, N)

tagged_reads = paired_reads(sampledpool, 250)
tagged_w_errs = mark_errors(tagged_reads, 0.001)

generate("tagged_reads.fastq", tagged_w_errs)

pool = Molecules("ecoli-ref.fasta", 5000)

cutter_a = fragment(40000)
tagger = tag(1000000)
cutter_b = fragment(700)
sampler = subsample(N) # Remember how to computed N previously.
mkreads = paired_reads(250)
adderr = mark_errors(0.001)

pool |> cutter_a |> tagger |> cutter_b |> sampler |> mkreads |> adderr |> generate("tagged_reads.fastq")

my_protocol = adderr ∘ mkreads ∘ sampler ∘ cutter_b ∘ tagger ∘ cutter_a

pool |> my_protocol |> generate("tagged_reads.fastq")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

