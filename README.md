# SequencingSim

*A package for simulating sequencing protocols and laboratory preparations.*

## Install

These instruction should get you going with this current WIP toolset:

1. Git clone this repo.
2. `cd` to the repo.
3. Use `julia --project -e'using Pkg; Pkg.instantiate()'`
4. `cd` to the scripts folder of the repo.
5. Use `julia --project -e'using Pkg; Pkg.instantiate()'`

Doing this should have made julia install the BioJulia and other packages this
project needs.

If this doesn't work... come yell at me...

## Available scripts

All scripts can be run using the julia command e.g.

`julia --project=scripts/ scripts/seqprep.jl --help`

The `--project` flag is important, set it to the directory containing these scripts.

To see the command line arguments for a script use the `--help` flag.

If any of the flags are unclear... yell at me to update the messages.

# 1. Create a workspace from a FASTA containing a genome

## a) In a julia session

```julia
using SequencingSim, BioSequences

open(FASTA.Reader, "mygenome.fasta") do rdr
    ws = workspace(rdr, false)
    save_workspace(ws, "myworkspace")
end
```

## b) Using `mkws_cmdline.jl`

```sh
julia mkws_cmdline.jl -I mygenome.fasta -O myworkspace
```


# 1. Copying up DNA sequences

The script `scripts/seqprep.jl` accepts an input FASTA file containing a genome
(like an a.lines.fasta of w2rap).

It makes many copies of each sequence. The number of copies of each sequence is
controlled by a mean and an SD, and the number of copies for a given genome is
drawn from a normal distribution.

You can also add a command line parameter to filter by size, only copying up
genome sequences of at least a given size.

It also needs an output file prefix.

# 2. Fragmentation of DNA sequences

The script `scripts/fragmenter_cmdline.jl` accepts an input FASTA file
containing sequences (typically prepared from script 1).

It fragments those input sequences and outputs a FASTA of the fragments.

It accepts a flag for the mean fragment length, and SD of the fragment length.

It also needs and output file prefix.

# 3. Sequencing of DNA fragments

The script `scripts/sequencer_cmdline.jl` accepts an input FASTA file of DNA
sequence fragments (typically prepared from script 2).

It sequences either end of the fragments and outputs a FASTQ file of the R1 reads,
and a FASTQ of the R2 reads.

You can set a mean and SD for the number of errors introduced per read. Which is then
drawn per read from a normal distribution (can't go lower than 0!).

The script requires an output prefix for output files.
It requires a fixed read length.
