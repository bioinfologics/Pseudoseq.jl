"""
    fabricate(file::String, cb...)

Fabricate the sequence(s) planned in a number chromosome blueprints.

The fabricated sequences will be written out to the FASTA formatted file `file`.

`cb...` should be provided as a series of blueprint, seed-sequence pairs.
E.g. `fabricate("mygenome.fasta", chr1plan => ch1seq, ch2plan => ch2seq)`.
"""
function fabricate(file::String, cb...)
    @info string("Opening FASTA file at: ", file)
    open(FASTA.Writer, file) do wtr
        fabricate(wtr, cb...)
    end
end

"""
    fabricate(fw::FASTA.Writer, cb...)

Fabricate the sequence(s) planned in a number chromosome blueprints.

The fabricated sequences will be written out to the FASTA formatted file `fw`.

`cb...` should be provided as a series of blueprint, seed-sequence pairs.
E.g. `fabricate("mygenome.fasta", chr1plan => ch1seq, ch2plan => ch2seq)`.
"""
function fabricate(fw::FASTA.Writer, cb...)
    chromno = 0
    for bp in cb
        sequences = fabricate(bp)
        chromno += 1
        for (i, seq) in enumerate(sequences)
            rec = FASTA.Record(string("Chromosome_", chromno, "_copy_", i), seq)
            write(fw, rec)
        end
    end
    @info string("Fabricated ", chromno, " chromosomes and written to file")
end

"""
    fabricate(cb::ChromosomeBlueprint)

Fabricate the sequence(s) planned in the chromosome blueprint.

A random DNA sequence will be generated to use as a seed sequence.

A sequence will be built for each chromosome copy in the blueprint.
"""
function fabricate(cb::ChromosomeBlueprint)
    seed = BioSequence{DNAAlphabet{2}}(randdnaseq(length(cb)))
    return fabricate(cb, seed)
end

"""
    fabricate(cb::ChromosomeBlueprint, seed::BioSequence{DNAAlphabet{2}})

Fabricate a DNA sequence by applying the planned features in a chromosome blueprint,
to some initial starting seed sequence.

A sequence will be built for each chromosome copy in the blueprint.
"""
function fabricate(cb::ChromosomeBlueprint, seed::BioSequence{DNAAlphabet{2}})
    apply_repeats!(seed, cb)
    copies = [copy(seed) for _ in 1:ncopies(cb)]
    apply_hets!(copies, cb)
    return copies
end
fabricate(p::Pair{ChromosomeBlueprint, BioSequence{DNAAlphabet{2}}}) = fabricate(p.first, p.second)
fabricate(p::Pair{BioSequence{DNAAlphabet{2}}, ChromosomeBlueprint}) = fabricate(p.second, p.first)

function apply_repeats!(seq::BioSequence{DNAAlphabet{2}}, cb::ChromosomeBlueprint)
    for repeat in repetitions_unsafe(cb)
        copyto!(seq, repeat.destination, seq, repeat.origin, repeat.size)
    end
    return seq
end

function apply_hets!(seqs::Vector{BioSequence{DNAAlphabet{2}}}, cb::ChromosomeBlueprint)
    for het in hets_unsafe(cb)
        @assert length(seqs) == length(het.alleles)
        p = het.position
        for (i, a) in enumerate(het.alleles)
            seqs[i][p] = a
        end
    end
end