using FASTX, BioSequences, Pseudoseq

refseq = open(FASTA.Reader, "yeast-chr1.fasta") do rdr
    FASTA.sequence(LongSequence{DNAAlphabet{2}}, read(rdr))
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

thet = plan_het(tetploid, .03, [1, 2], [3, 4])

thet = plan_het(thet, .01, [1, 3], [2, 4])

fabricate("build-a-yeast-tet.fasta", thet => refseq)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

