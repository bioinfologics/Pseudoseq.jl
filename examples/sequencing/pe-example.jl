# # Example: paired-end sequencing
#
# This is an example generated from this source
# file: [`pe-example.jl`](@__REPO_ROOT_URL__examples/sequencing/pe-example.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`pe-example.ipynb`](@__NBVIEWER_ROOT_URL__examples/sequencing/pe-example.ipynb)
#nb # jupyter notebook version. The corresponding online documentation page can
#nb # be found here: [`pe-example.html`](https://bioinfologics.github.io/Pseudoseq.jl/latest/examples/sequencing/pe-example),
# and the script can be found here: [`pe-example.jl`](./pe-example.jl) 


# For the simulation we are going to:
# 
# 1. Create a pool of 5000 copies of a reference genome.
# 2. Fragment the DNA molecules in the pool, to an average length of 700bp.
# 3. Subsample the molecules in the pool to achieve approximatly 50x coverage.
# 4. Create a set of 250bp paired-end reads.
# 5. Apply errors to the paired-end reads at a rate of 0.001 (.1%).
# 6. Generate an output FASTQ file.

using Pseudoseq.Sequencing

# ## Using the [`sequence`](@ref) method
# 
# First, let's see how we do this with the `sequence` method.
# The first two parameters we give to the function will be the input genome we
# want to sequence, and the destination FASTQ file for output reads.
# Here we are setting:
# - The number of genome copies in the molecule pool to 5000.
# - The average fragment size to 700bp.
# - The sampling coverage to 50x.
# - The read length to 250bp.
# - The per base read error rate to 0.001.
# - The fact we want paired-ends of fragments to be read (`paired`) to true.

sequence("ecoli-ref.fasta", "pe-reads.fastq"; ng = 5000, flen = 700, cov = 50, paired = true, rdlen = 250, err = 0.001)

# ## Using the `Molecules` tranformation methods.
# 
# Here's how to achieve the same thing, using the Pseudoseq API. It is nessecery to
# use the API if you want to do something that is not anticipated by the available
# functionality of the `sequence` method: the cost of conveinience is fewer options.
# 
# Starting with a FASTA formatted file containing the genome we want to sequence,
# we create a pool with 5000 copies of the genome.

pool = Molecules("ecoli-ref.fasta", 5000)

# Next we use the fragment function to make a pool of shorter DNA molecules.

cutpool = fragment(pool, 700)

# We need to determine the number of molecules to sample, and subsample the pool:

genome_size = 4639675
expected_coverage = 50
read_length = 250

N = needed_sample_size(expected_coverage, genome_size, read_length)
N = div(N, 2) # Divide by 2 as we're doing paired end sequencing.

sampledpool = subsample(cutpool, N)

# We now want to create a set of paired-end reads. We want our reads to be 250bp
# in length.

pe_reads = paired_reads(sampledpool, 250)

# Now we have some reads, we should mark positions in the reads that are destined
# to be errors in the output FASTQ.
# We will construct a `FixedProbSubstitutions` function with a per base error
# probability of 0.001 and pass it to the `edit_substitutions` method. 

f = FixedProbSubstitutions(0.001)
pe_w_errs = edit_substitutions(f, pe_reads)

# Now we have some paired end reads and have marked some positions as errors, we
# can generate FASTQ files.

generate("pe-reads.fastq", pe_w_errs)

# ## Constructing a pipeline of `Processors`.
#
# As a convenience, some users may prefer to use pipelines of `Processors`
# These behave like curried versions of the `Molecules` transformation methods.
# First let's define our starting `Molecules` pool:

pool = Molecules("ecoli-ref.fasta", 5000)

# To make a Processor, use a `Molecules` transformation method, but do not
# provide a `Molecules` value as a first argument. So let's make Processors for
# each step of our paired end sequencing pipeline.

cutter = fragment(700)
sampler = subsample(N) # Remember how to computed N previously.
mkreads = paired_reads(250)
adderr = mark_errors(0.001)

# Next we can construct the pipeline using standard julia function pipelining syntax:

pool |> cutter |> sampler |> mkreads |> adderr |> generate("pe-reads.fastq")

# You can also compose the processors together into one whole function.
# Typing \circ in the julia repl and then hitting tab gives you the circular
# composition symbol. Note how pipelining above progresses from left to right,
# but composition is right to left in order. 

my_protocol = adderr ∘ mkreads ∘ sampler ∘ cutter

pool |> my_protocol |> generate("pe-reads.fastq")

