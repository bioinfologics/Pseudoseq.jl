# 3. Subsampling the molecules

So we now have a universe of smaller DNA molecules. DNA sequencing experiments
are sampling processes, not every one of these fragmented molecules in the
universe will be read by the sequencer. So next we have to subsample the
universe.

This is done with the sample command:

```julia
subsampled_universe = sample(cut_universe, 463967)
```

Where the second numeric argument (463967) is the number of molecules you want
to sample.

## Determining the number of fragments to sample

But how did I get to the number 463967?

Prior to sampling, there is a universe of 33139992 DNA fragments, of
on average 700bp in length. There is variation - some fragments are bigger,
some are smaller.

We can estimate how often a base position in the genome is represented in this
universe of ~700bp fragments with a simple formula:

```math
C = \frac{L * N}{G}
```

Where `C` is the expected coverage, `L` is the average length of the fragments,
and `N` is the number of fragments, and `G` is the number of bases in the genome.

Given that `L = 700`, `N = 33139992`, and `G = 4639675`:

```math
C = \frac{700 * 33139992}{4639675}
```

The answer is ~5000, so each base in the genome is represented 5000 times 
(give or take) in the universe. This makes sense, because if you recall
how we created our universe.....

```julia
universe = makeuniverse("mygenome.fasta", 5000)
```

.....we started with 5000 copies of the genome!

Ok, so say I wanted to subsample the universe such that each base position in
the genome is represented ~50 times. That is, I want to achieve 50x coverage,
you can determine the number of fragments to subsample from the universe, by
reversing the formula:

```math
N = \frac{C * G}{L}
```

```math
N = \frac{50 * 4639675}{700}
```

The answer is 331405, so if I sample 331405 700bp fragments, from my universe
of 33139992 700bp fragments. I can expect that subsample to represent each base
in the genome approximately 50 times.

But wait! It's not quite a simple as that! The answer of 331405 would be correct
if we were to just pull those 700bp fragments out of the universe, and then
sequence a complete read of each 700bp fragment in the subsample in its entirety.

However, that's not what happens, even in long read sequencers. Rather some
sub-stretch of the fragment is read by the sequencer, and in the case of paired
end sequencing, a two sub-stretches of each fragment are read by the sequencer:
one from each end of the fragment:

```
Read 1
------->
TGAATAGCCCAAAGGTTTCGACA
|||||||||||||||||||||||
ACTTATCGGGTTTCCAAAGCTGT
              <--------
                 Read 2
```

And some of the fragment in the middle does not get read by the sequencer.

So how did I get the answer 463967?

First, I know that I want to simulate paired end sequencing, and that I want
the length of each read to be 250bp.

So I'm going to be sampling fragments of 700bp, from the universe, and reading
250bp in from either end.

I can estimate how many 700bp fragments I need to sample to achieve 50x coverage
by using the formula with a revised `L`, setting it to the read length 250bp,
instead of the fragment length of 700bp:

```math
N = \frac{50 * 4639675}{250}
```

This gives me an `N` of 927935, which I must divide by 2 (I'm going to be reading
two ends of each 700bp DNA fragment I sample), which gives me 463967.

Pseudoseq provides a helper function that assists in this type of calculation:

```julia
genome_size = 4639675
expected_coverage = 50
read_length = 250

N = needed_sample_size(expected_coverage, genome_size, readlength)

# Divide by 2 as we're doing paired end sequencing.
div(N, 2)
```

So after all that, we have a subsample of a universe of 700bp fragments.
The next step is to produce a set of reads from our universe sample.