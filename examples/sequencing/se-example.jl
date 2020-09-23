# # Example: long, single end reads
#
# This is an example generated from this source
# file: [`se-example.jl`](@__REPO_ROOT_URL__/examples/sequencing/se-example.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`se-example.ipynb`](@__NBVIEWER_ROOT_URL__/man/sequencing/examples/se-example.ipynb)
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

# Starting with a FASTA formatted file containing the genome we want to sequence,
# we create a pool with a copy of the genome.

m = Molecules("ecoli-ref.fasta")

# We want to first have 5000 full copies of the genome, so we will make an amplifier.

amp = amplify(5000)

# Next we want a fragmenter to make a pool of shorter DNA molecules, with an average
# length of 40,000bp.

frag = fragment(40000bp)

# Next we create a subsampler which will randomly sample molecules, to give us
# a desired expected coverage. We specify 40000bp as our predicted average read
# length.

ssmpl = subsample(50X, 40000bp)

# We then want a read-maker that will give us single end reads.
# Since we call `makereads` without providing a read length, the function will
# generate reads from the entire length of each molecule in the pool. We do this
# to emulate what Nanopore sequencing is supposed to do: It takes an entire DNA
# fragment, feeds it through an electrically charged pore, producing a read for
# the entire fragment.

readmaker = makereads()

# Once we have reads, we will mark positions in the reads that are incorrectly
# detected by the sequencer: errors.
# We will construct a `FixedProbSubstitutions` function with a per base error
# probability of 0.1 and pass it to the `make_substitutions` method. 
# This will make errors fall totally randomly over each read. 

errmaker = make_substitutions(FixedProbSubstitutions(0.1))

# Now we can push our molecules through a pipeline of these processors, and out
# to a FASTQ file:

m |> amp |> frag |> ssmpl |> readmaker |> errmaker |> generate("se-reads.fastq")

# You can also compose the processors together into one whole function.
# Typing \circ in the julia repl and then hitting tab gives you the circular
# composition symbol. Note how pipelining above progresses from left to right,
# but composition is right to left in order. 

my_protocol = errmaker ∘ readmaker ∘ ssmpl ∘ frag ∘ amp

m |> my_protocol |> generate("se-reads.fastq")