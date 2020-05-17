# # Example: long, single end reads
#
# This is an example generated from this source
# file: [`se-example.jl`](@__REPO_ROOT_URL__examples/sequencing/se-example.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`se-example.ipynb`](@__NBVIEWER_ROOT_URL__examples/sequencing/se-example.ipynb)
#nb # jupyter notebook version. The corresponding online documentation page can
#nb # be found here: [`se-example.html`](https://bioinfologics.github.io/Pseudoseq.jl/latest/examples/sequencing/se-example),
# and the script can be found here: [`se-example.jl`](./se-example.jl) 

# Let's see how you might simulate something like an Oxford Nanopore sequencing
# experiment.
#
# For the simulation we are going to:
# 
# 1. Create a pool of 5000 copies of a reference genome.
# 2. Fragment the DNA molecules in the pool, to an average length of 40,000bp.
# 3. Subsample the molecules in the pool to achieve approximatly 30x coverage.
# 4. Create a set of single-end reads, the enitre length of each molecule.
# 5. Apply errors to the reads at a rate of 0.10 (1 error every 10bp).
# 6. Generate an output FASTQ file.

using Pseudoseq.Sequencing

# ## Using the [`sequence`](@ref) method
# 
# First, let's see how we do this with the `sequence` method.
# The first two parameters we give to the function will be the input genome we
# want to sequence, and the destination FASTQ file for output reads.
# Here we are setting:
# - The number of genome copies in the molecule pool to 5000.
# - The average fragment size to 40000bp.
# - The sampling coverage to 30x.
# - The read length to `nothing`, which will make the sequencer read the whole length of any DNA fragment.
# - The per base read error rate to 0.1.
# - The fact we want paired-ends of fragments to be read (`paired`) to false.

sequence("ecoli-ref.fasta", "longreads.fastq"; ng = 5000, flen = 40000, cov = 30, rdlen = nothing, err = 0.1, paired = false)

# ## Using the Pseudoseq API
# 
# Here's how to achieve the same thing, using the Pseudoseq API. It is nessecery to
# use the API if you want to do something that is not anticipated by the available
# functionality of the `sequence` method: the cost of conveinience is fewer options.
# 
# Let's start with a pool of 5000 copies of a genome contained in a FASTA file:

pool = Molecules("ecoli-ref.fasta", 5000)

# Cut the pool of DNA into fragments of an average length of 40,000bp

cutpool = fragment(pool, 40000)

# Now we'll estimate the number of fragments we need to sample from the pool to
# achieve 30x coverage.

genome_size = 4639675
expected_coverage = 30
readlength = 40000

N = needed_sample_size(expected_coverage, genome_size, readlength)

sampledpool = subsample(cutpool, N)

# By using the [`make_reads`](@ref) function without specifying a read length,
# the function will generate reads from the entire length of each molecule in
# the pool. We do this to emulate what Nanopore sequencing is supposed to do:
# It takes an entire DNA fragment, feeds it through an electrically charged
# pore, producing a read for the entire fragment. 

se_reads = unpaired_reads(sampledpool, nothing)

# Long read sequencer have much higher error rates than short read sequencers
# so we use a error rate of 0.1 or roughly speaking an error every 10 bases.
# We will construct a `FixedProbSubstitutions` function with a per base error
# probability of 0.1 and pass it to the `edit_substitutions` method. 

f = FixedProbSubstitutions(0.1)
se_w_errs = edit_substitutions(f, se_reads)

# Finally produce the ouput FASTQ file.

generate("longreads.fastq", se_w_errs)

# ## Constructing a pipeline of `Processors`.
#
# As a convenience, some users may prefer to use pipelines of `Processors`
# These behave like curried versions of the `Molecules` transformation methods.
# First let's define our starting `Molecules` pool:

pool = Molecules("ecoli-ref.fasta", 5000)

# To make a Processor, use a `Molecules` transformation method, but do not
# provide a `Molecules` value as a first argument. So let's make Processors for
# each step of our single end sequencing pipeline.

cutter = fragment(40000)
sampler = subsample(N) # Remember how to computed N previously.
mkreads = unpaired_reads(nothing)
adderr = make_substitutions(FixedProbSubstitutions(0.1))

# Next we can construct the pipeline using standard julia function pipelining syntax:

pool |> cutter |> sampler |> mkreads |> adderr |> generate("se-reads.fastq")

# You can also compose the processors together into one whole function.
# Typing \circ in the julia repl and then hitting tab gives you the circular
# composition symbol. Note how pipelining above progresses from left to right,
# but composition is right to left in order. 

my_protocol = adderr ∘ mkreads ∘ sampler ∘ cutter

pool |> my_protocol |> generate("se-reads.fastq")