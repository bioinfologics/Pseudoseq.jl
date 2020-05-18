using Pseudoseq.Sequencing

sequence("ecoli-ref.fasta", "longreads.fastq"; ng = 5000, flen = 40000, cov = 30, rdlen = nothing, err = 0.1, paired = false)

pool = Molecules("ecoli-ref.fasta", 5000)

cutpool = fragment(pool, 40000)

genome_size = 4639675
expected_coverage = 30
readlength = 40000

N = needed_sample_size(expected_coverage, genome_size, readlength)

sampledpool = subsample(cutpool, N)

se_reads = unpaired_reads(sampledpool, nothing)

f = FixedProbSubstitutions(0.1)
se_w_errs = edit_substitutions(f, se_reads)

generate("longreads.fastq", se_w_errs)

pool = Molecules("ecoli-ref.fasta", 5000)

cutter = fragment(40000)
sampler = subsample(N) # Remember how to computed N previously.
mkreads = unpaired_reads(nothing)
adderr = make_substitutions(FixedProbSubstitutions(0.1))

pool |> cutter |> sampler |> mkreads |> adderr |> generate("se-reads.fastq")

my_protocol = adderr ∘ mkreads ∘ sampler ∘ cutter

pool |> my_protocol |> generate("se-reads.fastq")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

