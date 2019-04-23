using BioSequences, Pseudoseq

refseq = open(FASTA.Reader, "yeast-chr1.fasta") do rdr
    FASTA.sequence(BioSequence{DNAAlphabet{2}}, read(rdr))
end

reflen = length(refseq)

dploid = plan_chrom(reflen, 2)
dploidhet = plan_het(dploid, .01, 2)

# pair.

fabricate("build-a-yeast-di.fasta", dploidhet => refseq)

trploid = plan_chrom(reflen, 3)

trploidhet = plan_het(trploid, .01, 2)

# pair.

fabricate("build-a-yeast-tri.fasta", trploidhet => refseq)

tetploid = plan_chrom(reflen, 4)

thet = plan_het(tetploid, 6907, [1, 2], [3, 4])

import Random
Random.seed!(1234)
p = suggest_regions(tetploid, 1, 6907)
a1 = suggest_alleles(length(p), [1, 2], [3, 4])
a1 = suggest_alleles(length(p), [1, 1, 2, 2])
thet = plan_het(tetploid, 6907, [1, 2], [3, 4])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

