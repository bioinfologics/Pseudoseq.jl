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

A random motif is a motif that will change between calls to `make_puzzle`.
The reason is so as you can make multiple replicate puzzles with the same
properties, but with different specific sequences.

You can add a random motif to a `MotifStitcher` by using the `add_motif!` method,
and specifying only a length (in bp) for the random motif, for example:

```jlcon
ms = MotifStitcher
add_motif!(ms, 10_000)
```
