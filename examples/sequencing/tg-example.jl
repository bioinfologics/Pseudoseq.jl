# # Example: tagged paired-end reads
#
# This is an example generated from this source
# file: [`tg-example.jl`](@__REPO_ROOT_URL__examples/sequencing/tg-example.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`tg-example.ipynb`](@__NBVIEWER_ROOT_URL__examples/sequencing/tg-example.ipynb)
#nb # jupyter notebook version. The corresponding online documentation page can
#nb # be found here: [`tg-example.html`](https://bioinfologics.github.io/Pseudoseq.jl/latest/examples/sequencing/tg-example),
# and the script can be found here: [`tg-example.jl`](./tg-example.jl)

# Let's see how we might simulate something like an 10x sequencing experiment.
#
# For this simulation script we will:
#
# 1. Create a pool of 5000 copies of a reference genome.
# 2. Fragment the DNA molecules in the pool, to an average length of 40,000bp.
# 3. Tag the long molecules in the pool randomly with a set of 1,000,000 tags.
# 4. Fragment the molecules in the pool to an average length of 700bp.
# 5. Subsample the molecules in the pool to achieve approximatly 50x coverage.
# 6. Create a set of 250bp paired-end reads.
# 7. Apply errors to the paired-end reads at a rate of 0.001 (.1%).
# 8. Generate an output FASTQ file.

using Pseudoseq.Sequencing

# ## Using the [`sequence`](@ref) method
# 
# First, let's see how we do this with the `sequence` method.
# The first two parameters we give to the function will be the input genome we
# want to sequence, and the destination FASTQ file for output reads.
# Here we are setting:
# - The number of genome copies in the molecule pool to 5,000.
# - The number of possible tags to one million.
# - The average fragment size, prior to tagging, to 40,000bp.
# - The average fragment size after tagging, to 700bp.
# - The sampling coverage to 50x.
# - The read length to 250bp.
# - The per base read error rate to 0.001.
# - The fact we want paired-ends of fragments to be read (`paired`) to true, which is the default.

sequence("ecoli-ref.fasta", "tagged_reads.fastq"; ng = 5000, tusize = 1000000, taggedflen = 40000, flen = 700, cov = 50, rdlen = 250, err = 0.1)

# ## Using the Pseudoseq API
# 
# Here's how to achieve the same thing, using the Pseudoseq API. It is nessecery to
# use the API if you want to do something that is not anticipated by the available
# functionality of the `sequence` method: the cost of conveinience is fewer options.
#
# Let's start with a pool of 5000 copies of a genome contained in a FASTA file:

dnapool = Molecules("ecoli-ref.fasta", 5000)

# Now let's cut up the molecules to an average length of 40,000bp

cutpool = fragment(dnapool, 40000)

# Ok, now we will tag these large fragments randomly.
# Once you tag a fragment in a universe, any other fragments that
# are derived from that tagged fragment will inherit the same tag.

taggedpool = tag(cutpool, 1000000)

# Here I'm going to use a pool of 1,000,000 distinct tags.
# Which fragment gets a certain tag is random.
# The size of the tag pool, and the number of fragments in your universe will
# determine how likely it is that any two fragments get the same tag.
# Now we'll fragment the pool again 

taggedcutpool = fragment(taggedpool, 700)

# Subsample the pool of tagged molecules.

genome_size = 4639675
expected_coverage = 50
read_length = 250

N = needed_sample_size(expected_coverage, genome_size, read_length)
N = div(N, 2) # Divide by 2 as we're doing paired end sequencing.

sampledpool = subsample(taggedcutpool, N)

# Now let's make some 250bp tagged paired reads and generate some erroneous
# positions.
# We will construct a `FixedProbSubstitutions` function with a per base error
# probability of 0.001 and pass it to the `edit_substitutions` method. 

tagged_reads = paired_reads(sampledpool, 250)
f = FixedProbSubstitutions(0.001)
tagged_w_errs = edit_substitutions(f, tagged_reads)

# Output to FASTQ:

generate("tagged_reads.fastq", tagged_w_errs)

# ## Constructing a pipeline of `Processors`.
#
# As a convenience, some users may prefer to use pipelines of `Processors`
# These behave like curried versions of the `Molecules` transformation methods.
# First let's define our starting `Molecules` pool:

pool = Molecules("ecoli-ref.fasta", 5000)

# To make a Processor, use a `Molecules` transformation method, but do not
# provide a `Molecules` value as a first argument. So let's make Processors for
# each step of our paired-end, tagged-read, sequencing pipeline.

cutter_a = fragment(40000)
tagger = tag(1000000)
cutter_b = fragment(700)
sampler = subsample(N) # Remember how to computed N previously.
mkreads = paired_reads(250)
adderr = mark_errors(0.001)

# Next we can construct the pipeline using standard julia function pipelining syntax:

pool |> cutter_a |> tagger |> cutter_b |> sampler |> mkreads |> adderr |> generate("tagged_reads.fastq")

# You can also compose the processors together into one whole function.
# Typing \circ in the julia repl and then hitting tab gives you the circular
# composition symbol. Note how pipelining above progresses from left to right,
# but composition is right to left in order. 

my_protocol = adderr ∘ mkreads ∘ sampler ∘ cutter_b ∘ tagger ∘ cutter_a

pool |> my_protocol |> generate("tagged_reads.fastq")