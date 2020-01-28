# PuzzleMaker: _Core concepts & basic workflow_

`Pseudoseq` allows you to create puzzles; planned genomes and chromosome that
have a certain set of features and peculiarities of interest.

The purpose for creating such puzzle genomes is not to recreate biology perfectly.
The purpose is to create problems you understand fully (where the repeated
content is, which positions are heterozygous and so on).

Using such genomes can help you both understand and develop an intuition of what
current genome assembly tools are doing, and also to help design assembly
tools, and perhaps even plan sequencing experiments and form hypotheses.

`PuzzleMaker` is the `Pseudoseq` submodule that contains the functionality for
doing this. This manual includes several examples showing how to create genomes
with certain characteristics.
But the core workflow, and important concepts are explained below.

## The MotifStitcher

The `MotifStitcher` is the type that is central to `PuzzleMaker`.
It is by creating and interacting with a `MotifStitcher`, you create the plan
for your puzzle in 3 simple steps. First you define motifs, then you decide on
the order of motifs, and then you generate one or more instances of the puzzle
by calling `make_puzzle` on your `MotifStitcher`.

Let's take a look at step 1.

### 1. Define Motifs

The `MotifSticher` allows you to programmatically create puzzle haplotype
sequences by first defining a set of motifs, motifs are shorter chunks of user
specified or randomly generated sequence, that will be "stitched" together to
form the final haplotype sequences.

There are three kinds of motif that can be added to the MotifStitcher.

1. A random motif.
2. A fixed motif.
3. A sibling motif.

#### Random motifs

A random motif is a motif that will change between calls to `make_puzzle`.
The reason is so as you can trivially make multiple replicate puzzles with the
same properties, but with different specific sequences.

You can add a random motif to a `MotifStitcher` by using the `add_motif!` method,
and specifying only a length (in bp) for the random motif, for example:

```@setup pmkr
using Pseudoseq.PuzzleMaker
using BioSequences
```

```@repl pmkr
ms = MotifStitcher()
add_motif!(ms, 10_000)
```

By default, a random motif is constructed using a sampler that gives equal weighting
to the four nucleotides. If you wanted more control over some of the nucleotide
biases. You can construct a `RandomMotif` yourself, and provide it with your own
nucleotide sampler:

```@repl pmkr
ms = MotifStitcher()
smp = SamplerWeighted(dna"ACGT", [0.2, 0.3, 0.3,]) # Sampler biased toward GC.
add_motif!(ms, RandomMotif(10_000, smp)) # RandomMotif of 10_000 bp in length, using custom sampler.
```

!!! note
    Only `SamplerWeighted{DNA}` types are accepted.

#### Fixed motifs

Unlike a random motif, a fixed motif has its sequence defined and constant over
multiple calls of `make_puzzle`.

This is useful for situations where your puzzles must always include a certain
DNA sequence. Perhaps one of biological interest or known to confuse a heuristic.

You add a fixed motif to a `MotifStitcher` simply by passing it a DNA sequence:

```@repl pmkr
ms = MotifStitcher()
add_motif!(ms, dna"ATCGATCG")
```

#### Sibling motifs

A sibling motif is a motif that is randomly generated for each call of
`make_puzzle` just like random motifs. However, unlike random motifs, a sibling
motif is defined in terms of another motif already defined in the `MotifStitcher`.

To define a sibling motif, you specify an already existing motif. That motif's
sequence forms the base sequence of the new sibling motif. To define a sibling
motif you also need to provide a value that specifies the proportion of bases in
the new sibling motif's sequence, that should differ in their nucleic acid from
the base motifs sequence.

So sibling motifs then make it simple to define motifs that have a certain
level of sequence similarity / homology / shared ancestry, with another motif.
Creating portions of a simulated diploid genome might be one practical application
of sibling motifs.

You add a sibling motif to the `MotifStitcher` by providing the `add_motif!` method
with an `Pair{Int,Float64}`. where the integer is the ID of the chosen base motif
already defined in the `MotifStitcher`, and the floating point number specifies
the proportion of differing bases in the new motif's sequence:

```@repl pmkr
ms = MotifStitcher()
add_motif!(ms, 10_000) # A random first 10,000bp motif. Has ID = 1.
add_motif!(ms, 1 => 0.01) # Add a sibling motif that will have ~10 bases which differ from motif #1. Will have ID = 2.
```

### 2. Specify haplotypes

Once you have a set of motifs defined, you can build a set of haplotypes by
specifying sequences of motifs.

!!! note
    You specify the motifs using a vector of ID numbers. If you use a negative ID
    number -N then it means the reverse complement of the sequence of motif N.

!!! note
    A motif's ID can be repeated in such a vector any number of times, so you can
    create repeat structures in a haplotype. 

You add a haplotype by using the `add_motif_arrangement!` method with a
`MotifStitcher` and a vector of motif IDs.

For example, let's make a sequence that would form a hair-pin like structure,
with repeats when turned into a DeBruijn graph.

```@repl pmkr
ms = MotifStitcher()
add_motifs!(ms, 10000, 600, 10000, 600, 10000, 10000, 10000) # Use add_motifs! to add multiple random motifs at once.
add_motif_arrangement!(ms, [1, 2, 3, 4, 5, -4, -6, -2, 7])
```

### 3. Generate haplotype sequences

With the motifs defined and the haplotypes defined you can now generate sequences!

Simply call the `generate` method on the `MotifStitcher` to get a vector of
haplotype sequences. Repeatedly call `generate` to get independently generated
sequences from the same specification:

```@repl pmkr
generate(ms)
generate(ms)
```

If you provide a filename, the haplotypes will be written to file in FASTA
format, instead of returned as a value:

```@repl
generate(ms, "myhaplos.fasta")
```

That's all there is to it. Now you can try a simulated sequencing experiment on
your haplotypes. 