# # Example: tagged paired-end reads
#
# This is an example generated from this source
# file: [`tg-example.jl`](@__REPO_ROOT_URL__/examples/sequencing/tg-example.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`tg-example.ipynb`](@__NBVIEWER_ROOT_URL__/man/sequencing/examples/tg-example.ipynb)
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

# Let's start with a pool of 5000 copies of a genome contained in a FASTA file:

m = Molecules("ecoli-ref.fasta")

# We want to first have 5000 full copies of the genome, so we will make an amplifier.

amp = amplify(5000)

# Next we want a fragmenter to make a pool of shorter DNA molecules, with an average
# length of 40,000bp.

firstfrag = fragment(40000bp)

# Ok, now we want to assign random tags to these large fragments.
# Once you tag a molecule any other molecules that are derived from that tagged
# molecule will inherit the same tag.
# Our tagger will have a pool of 1,000,000 distinct tags.
# Which fragment gets a certain tag is random.
# The size of the tag pool, and the number of fragments in your universe will
# determine how likely it is that any two fragments get the same tag.

tagger = tag(1000000)

# We'll want to fragment the molecules again to even shorter lengths. 

secondfrag = fragment(700bp)

# We will want a filter to make sure no extreme sized molecules make their way through.

size_filter = select(x -> 900 >= length(x) >= 450)

# Next we create a subsampler which will randomly sample molecules, to give us
# a desired expected coverage. We will want our reads to be 250bp in length.

ssmpl = subsample(50X, 2 * 250bp)

# We then want a read-maker that will give us paired end reads. We want our
# reads to be 250bp in length.

readmaker = makereads(2 * 250bp)

# Once we have reads, we will mark positions in the reads that are incorrectly
# detected by the sequencer: errors.
# We will construct a `FixedProbSubstitutions` function with a per base error
# probability of 0.001 and pass it to the `make_substitutions` method. 
# This will make errors fall totally randomly over each read. 

errmaker = make_substitutions(FixedProbSubstitutions(0.001))

# Now we can push our molecules through a pipeline of these processors, and out
# to a FASTQ file:

m |> amp |> firstfrag |> tagger |> secondfrag |> size_filter |> ssmpl |> readmaker |> errmaker |> generate("tagged_reads.fastq")

# You can also compose the processors together into one whole function.
# Typing \circ in the julia repl and then hitting tab gives you the circular
# composition symbol. Note how pipelining above progresses from left to right,
# but composition is right to left in order. 

my_protocol = errmaker ∘ readmaker ∘ ssmpl ∘ size_filter ∘ secondfrag ∘ tagger ∘ firstfrag ∘ amp

m |> my_protocol |> generate("tagged_reads.fastq")