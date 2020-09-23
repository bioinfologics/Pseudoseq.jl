# # Example: paired-end sequencing
#
# This is an example generated from this source
# file: [`pe-example.jl`](@__REPO_ROOT_URL__/examples/sequencing/pe-example.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`pe-example.ipynb`](@__NBVIEWER_ROOT_URL__/man/sequencing/examples/pe-example.ipynb)
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

# Starting with a FASTA formatted file containing the genome we want to sequence,
# we create a pool with a copy of the genome.

m = Molecules("ecoli-ref.fasta")

# We want to first have 5000 full copies of the genome, so we will make an amplifier.

amp = amplify(5000)

# Next we want a fragmenter to make a pool of shorter DNA molecules, with an average
# length of 700bp.

frag = fragment(700bp)

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

m |> amp |> frag |> size_filter |> ssmpl |> readmaker |> errmaker |> generate("pe-reads.fastq")

# You can also compose the processors together into one whole function.
# Typing \circ in the julia repl and then hitting tab gives you the circular
# composition symbol. Note how pipelining above progresses from left to right,
# but composition is right to left in order. 

my_protocol = errmaker ∘ readmaker ∘ ssmpl ∘ size_filter ∘ frag ∘ amp

m |> my_protocol |> generate("pe-reads.fastq")