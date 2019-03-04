# # Example: paired-end sequencing
#
# This is an example generated from this source
# file: [`pe-example.jl`](@__REPO_ROOT_URL__examples/pe-example.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`pe-example.ipynb`](@__NBVIEWER_ROOT_URL__generated/pe-example.ipynb)
#nb # jupyter notebook version. The corresponding online documentation page can
#nb # be found here: [`pe-example.html`](https://bioinfologics.github.io/dev/generated/pe-example.html),
# and the script can be found here: [`pe-example.jl`](./pe-example.jl) 


# For the simulation we are going to:
# 
# 1. Create a pool of 5000 copies of a reference genome.
# 2. Fragment the DNA molecules in the pool, to an average length of 700bp.
# 3. Subsample the molecules in the pool to achieve approximatly 50x coverage.
# 4. Create a set of 250bp paired-end reads.
# 5. Apply errors to the paired-end reads at a rate of 0.001 (.1%).
# 6. Generate an output FASTQ file.

using Pseudoseq

# Starting with a FASTA formatted file containing the genome we want to sequence,
# we create a pool with 5000 copies of the genome.

pool = makepool("../ecoli-ref.fasta", 5000)

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

pe_reads = make_reads(PairedEnd, sampledpool, 250)

# Now we have some reads, we should mark positions in the reads that are destined
# to be errors in the output FASTQ.

pe_w_errs = mark_errors(pe_reads, 0.001)

# Now we have some paired end reads and have marked some positions as errors, we
# can generate FASTQ files.

generate("pe-reads.fastq", pe_w_errs)

# That's all there is to it!