# 1. Creating the universe

All DNA sequencing experiments begin with a sample of tissue or cells. Some hair,
some blood and so forth.

Such a sample undergoes a DNA extraction preparation in the laboratory, after which
the genome exists as number of DNA molecules, suspended in a solution.

Many cells are typically used as raw input material for DNA extraction, and so the
extracted DNA material contains a great many copies of the genome.

In `Pseudoseq`, we call this extracted genetic material (or rather the abstraction
of it) *the universe*. So named, as it is the totality of all copies of the
genome that exist in your simulation, from which all subsequent library prep will
be done, and from which all reads will be sequenced.

`Pseudoseq` abstracts DNA sequencing experiments as sampling processes,
because this is what they are from a statistical point of view: Just as a quadrat
placed at random on the forest floor provides a small sample of it's species
composition, so it is that a sequencing read provides a small sample of the
composition of motifs present in a genome.
This sample is the universe from which all samples will be drawn.

Starting with a FASTA formatted file containing the genome that you want to simulate
sequencing for, we create such a universe with `Pseudoseq` as follows:

```julia
universe = makeuniverse("mygenome.fasta", 5000)
```

Where the second argument is the number of copies of the genome you want to
exist in your universe. The example above would create a universe of 5000 copies
of the genome in "mygenome.fasta".

You can think of your universe at this stage as containing many copies of the
same genome sequence:

```
   1. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
   2. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
   3. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
                                    ...
4999. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
5000. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
```