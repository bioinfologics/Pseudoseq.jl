# Build-a-Genome: _Core concepts & basic workflow_

`Pseudoseq` allows you to plan and build genomes and chromosomes that have
a certain set of features and peculiarities. The purpose for doing this is
not to recreate biology perfectly. The purpose is to create genomes you understand
fully (where the repeated content is, which positions are heterozygous and so on).

Using such genomes can help you both understand and develop an intuition of what
current genome assembly tools are doing, and also to help design assembly
tools, and perhaps even plan sequencing experiments and form hypotheses.

This manual includes several examples showing how to plan genomes with certain
characteristics. But the core workflow, and important concepts are explained below,
in 3 steps.


## 1. Creating chromosome blueprints

Chromosome blueprints are the backbone of simulating genomes with Pseudoseq.

Chromosome blueprints determine the nature of one chromosome in a genome.

You can think of chromosome blueprints in Pseudoseq as: _a collection of
operations which, when applied to some seed sequence, result in a set of N
homologous sequences_.

!!! note
    Chromosome blueprints are _immutable_: Adding an operation to a chromosome
    blueprint creates a new blueprint, which is a copy of the old blueprint with
    the addition of the new operation.

!!! note
    Ideally, for any chromosome blueprint you construct with Pseudoseq, for each
    operation is possesses, it must be possible (at least in principle) to be able to
    intuit what the effect of that operation will be on:
    
    1. The Khmer distribution produced by sequencing reads of the fabricated chromosome.
    2. The structure of a sequence graph produced by sequencing reads of the fabricated chromosome.
    
    Therefore, we have made the design decision that _no two operations in a
    chromosome blueprint may affect the same position(s) of the genome in a
    conflicting manner._ To meet this requirement, certain operations "consume"
    a region of the chromosome planned in a blueprint. If a region is consumed,
    another operation that would also affect that region cannot be added to the
    blueprint.

Depending on the genome, any given chromosome may be present in multiple copies.
Diploids, for example have two copies of every chromosome.

The first step in simulating any artificial genome is to create one or more
blank chromosome blueprints. The [`plan_chrom`](@ref) function is used for this.

For example, this:

```@repl
using Pseudoseq
c = plan_chrom(100, 2)
```

will create a blank blueprint for 2 copies of a chromosome of 100bp length.
From the output you can see that it prints for you the number of chromosome copies,
the length of the chromosomes, and a list of available, _unconsumed_ regions of
the chromosome (see note above).

## 2. Adding planned features to a chromosome blueprint

After creating a fresh chromosome blueprint, no plans (operations) have been
added yet.

If you were to [`fabricate`](@ref) this blank blueprint, you would get `N`
identical DNA sequences as output, where `N` is the number of copies the
blueprint was planning.

Once you have one or more chromosome blueprints, you can add features to them.

This is done with a series of consistently named `plan_*` functions.

!!! note
    Remember; blueprints are immutable, so every time one of these `plan_*`
    functions is used to add a feature to a chromosome blueprint, a new chromosome
    blueprint is created.

### Repetitions

A repetition is a segment of a sequence, that has the exact same motif, as 
another segment of the sequence.

In Pseudoseq, to plan a repetition, you specify a region the repetition will copy,
sometimes called the `from` region. You also specify a region where the motif in
the `from` region will be replicated, called the `to` region.
So you might find it helpful to imagine a planned repetition as a kind of
copy-paste operation that occurs during [`fabricate`](@ref).

!!! note
    Repetitions consume the `to` region of the chromosome blueprint to which they
    are added.
    Repetitions _do not_ consume the `from` region, so other operations are free
    to affect the motif in the `from` region.
    Just remember that the repetition will replicate anything in the `from` region,
    including other features such as heterozygosity that occur in `from`.

Repetitions are added to a genome blueprint using the [`plan_repetition`](@ref) function.

### Heterozygosity

A heterozygosity describes a base position at which the copies of the chromosome
in the blueprint differ.

For a blueprint with 2 copies, both copies will differ at a given position.

For a triploid, at a heterozygous position, all 3 copies might differ from each
other. Alternatively, it is possible that 2 copies are the same, but they differ
from a 3rd copy. This applies for blueprints with higher copy numbers too: at
a heterozygous position some copies will differ from each other at that position,
but some of the copies might be the same.

!!! note
    Heterozygosity operations consume the position at which they are defined.

Heterozygous positions are planned using the [`plan_het`](@ref) function.

For example, below will make it so as about 50% of the two chromosome copies 
are heterozygous:

```@repl
using Pseudoseq
c = plan_chrom(100, 2)
chet = plan_het(c, .50, 2)
```

The above use of [`plan_het`](@ref) allows the function to choose which sites in the
blueprint are heterozygous, and how to allocate the bases at the heterozygous
sites. We only instruct the function that 50% of the genome should be heterozygous,
and that there should be 2 possible bases at each position (the only option really:
the blueprint plans for a diploid - 2 chromosome copies).

You can also take more fine-grained control:

```@repl het
using Pseudoseq, BioSequences
c = plan_chrom(100, 2)
chet = plan_het(c, 20, [DNA_T, DNA_A])
```

Here I planned a heterozygous position, at one site in the chromosome (20), and
I set the state of the first copy to `DNA_T`, and the state of the second copy to
`DNA_A`.

So as you can see, [`plan_het`](@ref) is very flexible. Check it's API documentation,
and the "Build-A-Yeast" example walkthrough to see more examples of [`plan_het`](@ref)
use.

### Utility functions

There is a utility function [`suggest_regions`](@ref) available to help you plan
where to place features.

Say you wanted to see where you could place 3 5bp repetitions, you
could do the following:

```@repl het
r = suggest_regions(chet, 5, 6)
```

So now you have 6 5bp regions, every 2 regions defining a `from` and a `to` region
for a single repetition. 

You can provide `r` as an input to [`plan_repetition`](@ref).

```@repl het
crep = plan_repetition(chet, r)
```

Another utility function [`suggest_alleles`](@ref) is available to help you plan
nucleotide patterns at heterozygous sites.

Say you wanted to plan a pattern in which two of three chromosome copies had the
same base, and a third differed.

```@repl het
suggest_alleles(3, 2)
```

The above asks for an allele pattern for 3 copies of a chromosome, with 2 possible
alleles, which chromosome copy gets which allele is random.

You can also take more control and specify which chromosome copy gets which state:

```@repl het
suggest_alleles(chet, [1, 1, 2])
# Or alternatively put:
suggest_alleles(chet, [1, 2], [3])
```

In the above, which chromosomes get the same allele is user-specified, but which
bases those alleles are, is randomly determined.

You can get many suggestions at once:

```@repl het
suggest_alleles(5, [1, 1, 2])
```

See the API docs for [`suggest_alleles`](@ref) for more details on the arguments
permitted by the different methods. 

## 3. Fabricate a FASTA from chromosome blueprints

Once you have a set of chromosome blueprints with the features planned that you
desire, you can fabricate the sequences of these chromosomes.

To do that, use the [`fabricate`](@ref) method.

You can simply fabricate the sequences, for use in the interactive session:

```@repl het
fabricate(chet)
```

Or have the sequences output to a FASTA formatted file.

```@repl het
fabricate("mychrom.fasta", chet)
```

A randomly generated seed sequence is used to start the fabrication process
unless you provide one. See the API docs for [`fabricate`](@ref) for more details.

And that's all there is to building a chromosome with Pseudoseq.
To build a genome with more than one chromosome, simply build a set of chromosome
blueprints.

Now check out the examples to see how various genomes can be built with Pseudoseq.
