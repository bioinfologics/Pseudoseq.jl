# # Example: tagged paired-end reads
#
# This is an example generated from this source
# file: [`tg-example.jl`](@__REPO_ROOT_URL__examples/tg-example.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`tg-example.ipynb`](@__NBVIEWER_ROOT_URL__generated/tg-example.ipynb)
#nb # jupyter notebook version. The corresponding online documentation page can
#nb # be found here: [`tg-example.html`](https://bioinfologics.github.io/dev/generated/tg-example.html),
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

using Pseudoseq

dnapool = makepool("../ecoli-ref.fasta", 5000)
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

tagged_reads = make_reads(TaggedPairs, sampledpool, 250)
tagged_w_errs = mark_errors(tagged_reads, 0.001)

# Output to FASTQ:

generate("tagged_reads.fastq", tagged_w_errs)