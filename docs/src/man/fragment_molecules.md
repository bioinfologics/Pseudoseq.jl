# 2. Fragmenting the molecules

So we have created a perfect universe of 5000 genome molecules.

If DNA sequencing machines could start at one end of a molecule and sequence
all the way to the end with reasonable accuracy and throughput, all we would have
to do is select a full DNA molecule from our perfect universe, and run it through
a sequencing process to get our genome read. Actually, we'd run several of the
molecules from our perfect universe through the sequencer, and use those 
reads to form a consensus that eliminates any mistakes a sequencer machine might
make - no technology is completely error proof after all and all detection
equipment has a certain rate of error.

Sadly, we don't have such an amazing sequencer technology that allows the reading
of a chromosome in its entirety. Instead, shorter reads are taken of shorter
fragments of DNA (even long read technology uses fragments and produces reads,
much shorter than the full size of a chromosome molecule).

A common step in any DNA sequencing experiment, therefore, is to fragment or
shear the DNA molecules that are present in the sample.

So a common step in any sequencing experiment simulation constructed with
`Pseudoseq` is to fragment the DNA molecules in a universe into smaller ones.

This is achieved with the `fragment` method:

```julia
cut_universe = fragment(universe, 700)
```

The second value of 700 provided above is the desired expected length of the
fragments. The `fragment` method uses this value, and the length of each
molecule it is about to cut, to decide how many breakpoints to scatter across
said fragment.

You can could visualise this process like so:

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