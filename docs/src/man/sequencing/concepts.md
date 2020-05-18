# Sequencing: _Core concepts & basic workflow_

`Pseudoseq` abstracts DNA sequencing experiments as sampling processes,
because this is what they are from a statistical point of view.
Just as a quadrat placed at random on the forest floor provides a small sample
of it's species composition, so it is that a sequencing read provides a small
sample of the composition of the motifs present in a genome.

This manual includes several examples showing how to emulate various sequencing
experiments using different technologies. But the core workflow, and important
concepts are outlined below.

!!! tip
    The user can use Pseudoseq's API to script each stage of the flow outlined
    below themselves, or they can use the [`sequence`](@ref) function, which is
    the highest-level user facing function. Every example in this section of
    the manual, will show you how to use both options to achieve the same goal.

Anyway, the core sequencing workflow in Pseudoseq is as follows...

## 1. Create a pool of DNA molecules

In reality, all DNA sequencing experiments begin with a sample of tissue or cells.
DNA is extracted from the sample in the laboratory, after the genomes of the
cells exist as number of DNA molecules, suspended in a solution.

In `Pseudoseq`, such a collection of DNA molecules is called *a molecule pool*,
it is created with the [`Molecules`](@ref) constructor function.

In the beginning of a `Psueodseq` simulation script you create a pool that is
the totality of all copies of the genome that exist in your simulation.

!!! tip
    You can think of the pool at this stage as containing many copies of
    the same genome sequence:
    
    ```
       1. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
       2. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
       3. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
                                    ...
    4999. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
    5000. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
    ```


## 2. Process the DNA molecule pool

You then subject your starting DNA molecule pool to a series of transformations,
until the pool has the properties you want to emulate.

Here we will describe the different transformations you can apply to a DNA
molecule pool:


### Fragmenting the pool

In an ideal world, if DNA sequencing machines could start at one end of a
molecule and read the sequence all the way to the end with reasonable accuracy
and throughput, you could simulate that by selecting molecules from your pool,
and producing a read file. Give or take some inevitable errors (all detection
equipment has a rate of error), assembling a genome would be simple.

Sadly, we don't have such an amazing sequencer technology that allows the reading
of a chromosome in its entirety. Instead, shorter reads are taken of short
fragments of DNA.
Even long read technology uses fragments and produces reads, much shorter than
the full size of a chromosome molecule.

A common step in any DNA sequencing experiment, therefore, is to fragment or
shear the DNA molecules that are present in the pool.

This is achieved in `Pseudoseq` with the [`fragment`](@ref) function.

!!! tip
    You can visualise this process like so:
    
    From:
    ```
       1. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
       2. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
       3. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
                                    ...
    4999. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
    5000. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC
    ```
    
    To:
    ```
       1. CGGACTT GAATAGC CCAAA GGTTTCGACACGA TCACGAC ACATAAAT TGGCGGAC TTGAATAGC
       2. CGGA CTTGAAT AGCCCAAAG GTTTCGAC ACGATCACGACACAT AAATTGGCGGA CTTGA ATAGC
       3. CGGACTTGA ATAGCC CAAAGGT TTCGACACGAT CACGACACA TAAATT GGCGGACTT GAATAGC
                                    ...
    4999. CGGAC TTGAATA GCCCAAAGGTTT CGACACGA TCACGACACAT AAATTG GCGGACTTG AATAGC
    5000. CGGACTTGAA TAGCCCA AAGGTTTCGA CACGATCAC GACACA TAAATTGGCGG ACTTGAAT AGC
    ```


### Subsampling molecules from a pool

DNA sequencing experiments are sampling processes, not every one of the 
DNA molecules in a pool will be read by the sequencer.

Randomly subsampling (without replacement) DNA molecules of a pool is achieved
in `Pseudoseq` with the [`subsample`](@ref) function.


#### Determining the number of molecules to subsample

You can estimate how often a base position in the genome is represented in a
DNA molecule pool using the following formula.

```math
C = \frac{LN}{G}
```

Where ``C`` is the expected coverage, ``L`` is the average length of the fragments,
and ``N`` is the number of fragments, and ``G`` is the number of bases in the genome.

If you wanted to subsample a pool such that each base position in a genome is
represented ~50 times (i.e. to achieve 50x coverage), you can determine the
number of fragments to subsample from the universe, by reversing the formula:

```math
N = \frac{CG}{L}
```

For example, if you wanted to sequence 250bp paired-end reads, at 50x coverage,
with a genome size of 4639675bp:

```math
\frac{50 \times 4639675}{250} = 927935
```

Remembering that you get two reads from one DNA molecule with paired-end sequencing,
you know to subsample ``927935 / 2 = 463967`` DNA molecules from a pool.

!!! tip
    Pseudoseq provides a helper function that assists in this type of calculation:
    
    ```julia
    genome_size = 4639675
    expected_coverage = 50
    read_length = 250
    
    N = needed_sample_size(expected_coverage, genome_size, readlength)
    
    # Divide by 2 as we're doing paired end sequencing.
    div(N, 2)
    ```


### Tagging molecules in a pool

Some DNA sequencing technologies work by allowing short reads to contain longer
range information by tagging molecules in a pool. The basic idea being that if
two short reads come from the same large DNA molecule, they will have the same
tag.

If a DNA fragment in a universe is tagged, and then it is subsequently fragmented
during a [`fragment`](@ref) transform, then all the smaller fragments derived
from that long fragment will inherit that long fragment's tag. This allows
shorter fragments to possess longer range information in the form of these tags,
and this is the basis of 10x sequencing and similar technologies.

`Pseudoseq` lets you attach tags to molecules in a pool using the [`tag`](@ref)
function.


## 3. Generating reads

Next, you generate a set of reads from your transformed pool. 
`Pseudoseq` allows you to create paired-end sequencing reads, single-end
sequencing reads, and linked paired-end sequencing reads.

You create a set of reads using the [`paired_reads`](@ref) or
[`unpaired_reads`](@ref) functions.

At this point in the process, the reads will be perfect. Generating perfect
reads might be useful for some purposes, such as for education: illustrating 
assembly processes initially without the complication of errors to aid comprehension.
It might also be interesting to illustrate the action of a heuristic on data without
errors cf. data that does have errors. 

In reality, all sequencers have certain error profiles - certain genomic motifs
can cause different technologies to fail in different ways, and the probability
of errors can change over time, particularly if the sequencing process involves
chemical reactions, or if synchronisation is important. Some Illumina sequencers,
for example, are subject to what is known as the phase-problem.

So how do we create errors?

Pseudoseq currently supports errors in the form of single base pair substitutions.

The reads object contains a Dictionary of (`Integer` => `Vector{Substitution}`)
pairs, where the integer is the (1 based) index of the read, and the vector of
substitutions contains the substitutions that will be applied to the sequence of
the read at read-file generation.

The substitution type Pseudoseq uses is very simple can can be created using a
base position, nucleotide pair. For example `Substitution(51, DNA_A)` results
in a subsitution which would replace base 51 of a read with an 'A' nucleotide.

Whilst you could populate and edit this dictionary manually, doing it for every
read would be tedious.

You can generate substitutions programmatically though using the
[`edit_substitutions`](@ref) method, which allows you to pass an anonymous
function or "functor type", which edits the substitutions for each read.

Pseudoseq provides some built in functions for use with [`edit_substitutions`](@ref),
and more will be added over time. A list with links to their API documentation
is provided here.

1. [`FixedProbSubstitutions`](@ref)
2. [`ClearSubstitutions`](@ref)

Please see the API docs for [`edit_substitutions`](@ref) for details of the
requirements for any functions passed to [`edit_substitutions`](@ref).

If you have any questions about how to model a particular error profile, please
see the examples section of the documentation, and also don't hesitate to open a
GitHub Issue describing your particular scenario.

Finally, you use your set of reads to [`generate`](@ref) either an interleaved
FASTQ file, or two FASTQ files (one for R1 reads, and one for R2 reads).