# # Example: Build-a-Yeast
#
# This is an example generated from this source
# file: [`build-a-yeast.jl`](@__REPO_ROOT_URL__examples/build-a-genome/build-a-yeast.jl)
# You are seeing the
#md # online documentation version. The corresponding notebook can be found
#md # here: [`build-a-yeast.ipynb`](@__NBVIEWER_ROOT_URL__examples/build-a-genome/build-a-yeast.ipynb)
#nb # jupyter notebook version. The corresponding online documentation page can
#nb # be found here: [`build-a-yeast.html`](https://bioinfologics.github.io/Pseudoseq.jl/dev/examples/build-a-genome/build-a-yeast),
# and the script can be found here: [`build-a-yeast.jl`](./build-a-yeast.jl) 

# We are going to use Pseudoseq to create a simple set of fake genomes, based on
# chromosome 1 of the reference genome for yeast.
# 
# 1. Load the sequence of chromosome 1 of the yeast reference from FASTA file.
# 2. Create blank chromosome blueprints for different genome designs.
# 3. Edit the blueprints to achive a:
#    a. Diploid chromosome
#    b. Triploid chromosome
#    c. Tetraploid chromosome

using BioSequences, Pseudoseq

# ## Load the seed sequence
# 
# First load the sequence from FASTA file. This uses tools from BioSequences,
# which is a dependency of Pseudoseq and so the user should have it available,
# as julia's package manager would have installed any dependencies.
# The lines below open a fasta file, with a fasta reader, and load a single
# FASTA record, and get its sequence. 

refseq = open(FASTA.Reader, "yeast-chr1.fasta") do rdr
    FASTA.sequence(BioSequence{DNAAlphabet{2}}, read(rdr))
end

reflen = length(refseq)

# This results in a 230,218nt long sequence called `refseq`.
# 
# ## Designing a diploid
# 
# First lets make a plan for a diploid genome with 1% heterozygosity.
# Pseudoseq has a lot of flexibilty in how you can define the heterozygosity,
# but for ease here I'm going to let it decide where the sites should be, and
# what the alleles at each site should be, by calling [`plan_het`](@ref),
# providing just the proportion of sites I want to be heterozygous, and the
# number of alleles I want at each of these heterozygous sites (2).

dploid = plan_chrom(reflen, 2)
dploidhet = plan_het(dploid, .01, 2)

# Now fabricate the FASTA, with [`fabricate`](@ref), providing a filename, and then a
# series of `blueprint => seed` pairs. In this case we just give 1 blueprint-seed
# pair.

fabricate("build-a-yeast-di.fasta", dploidhet => refseq)

# Thats it for the diploid genome.
# 
# ## Designing a triploid
# 
# Now let's make a triploid genome. This should also be simple enough. It's much
# like the process for making the diploid genome, but specifying 3 chromosome
# copies instead of 2.

trploid = plan_chrom(reflen, 3)

# When I add some heterozygous sites, just like with the diploid case, I'm going
# to simply ask for a proportion of heterozygous sites (and let Pseudoseq decide
# where they should occur. Like last time also I'm going to make all the sites have 2
# alleles, so for every heterozygous site, two of the three chromosome copies
# willl have the same base, and one of the copies will be different. There are
# other ways to define some heterozygous sites for example all 3 chromosome copies
# could have different bases

trploidhet = plan_het(trploid, .01, 2)

# Now fabricate the FASTA, with `fabricate`, providing a filename, and then a
# series of `blueprint => seed` pairs. In this case we just give 1 blueprint-seed
# pair.

fabricate("build-a-yeast-tri.fasta", trploidhet => refseq)

# That's it for the triploid genome.
# 
# ## Designing a tetraploid genome
# 
# Lets make a tetraploid genome, that consists of two diploid genomes, an A
# genome and a B genome.
# First I'm going to create a chromosome blueprint with 4 copies:

tetploid = plan_chrom(reflen, 4)

# We want to plan the heterozygosity in this tetraploid such that the
# number of heterozygous sites within the A and B genomes is ~1%. But we want
# ~3% of sites to be heterozygous between the A and B genomes.
# 
# First we will say copies 1 + 2 in our blueprint are the A genome, and copies
# 3 + 4 are the B genome. Next consider heterozygous sites with 2 different alleles,
# A and B, 6 different patterns are possible:
# 
# | Copy 1 | Copy 2 | Copy 3 | Copy 4 |
# |:------:|:------:|:------:|:------:|
# | A      | A      | B      | B      |
# | B      | B      | A      | A      |
# | A      | B      | A      | B      |
# | B      | A      | B      | A      |
# | A      | B      | B      | A      |
# | B      | A      | A      | B      |
# 
# The first 2 patterns in the table constitute fixed mutations between 
# genomes A and B,
# but not polymorphism within genomes A and B. This is because Copy 1 and 3 differ,
# as do copies 2 and 4, yet copies 1 and 2 are the same, as are copies 3 and 4.
# To achieve a ~3% divergence between genomes A and B, let's plan some
# divergent mutations by grouping the copy numbers together when planning some
# heterozygous sites:

thet = plan_het(tetploid, .03, [1, 2], [3, 4])

# By using the groups `[1, 2]` and `[3, 4]` we ensure chromosome copy 1 and 2
# will have the same allele, as will chromosomes 3 and 4

# Next we'll just add some of the other heterozygosity patterns to achieve the
# 1% internal heterozygosity desired. We'll use a heterozygosity pattern with
# two states.

thet = plan_het(thet, .01, [1, 3], [2, 4])

fabricate("build-a-yeast-tet.fasta", thet => refseq)
