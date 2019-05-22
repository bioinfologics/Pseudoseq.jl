
struct MoleculePool
    genome::Vector{BioSequence{DNAAlphabet{2}}}
    seqviews::Views
end

"""
    makepool(gen::Vector{BioSequence{DNAAlphabet{2}}}, ng::Int = 1)

Create a pool of `ng` copies of a genome defined by the `gen` vector of sequences.
"""
function makepool(gen::Vector{BioSequence{DNAAlphabet{2}}}, ng::Int = 1)
    vs = Views()
    for (i, seq) in enumerate(gen)
        push!(vs, SequencingView(i, 1, length(seq), false))
    end
    if ng > 1
        vs = amplify(vs, ng)
    end
    return MoleculePool(gen, vs)
end

"""
    makepool(rdr::FASTA.Reader, ng::Int = 1)

Create a pool of `ng` copies of the genome read in from the `FASTA.Reader`.
"""
function makepool(rdr::FASTA.Reader, ng::Int = 1)
    gen = Vector{BioSequence{DNAAlphabet{2}}}()
    rec = eltype(rdr)()
    while !eof(rdr)
        read!(rdr, rec)
        push!(gen, FASTA.sequence(BioSequence{DNAAlphabet{2}}, rec))
    end
    return makepool(gen, ng)
end

"""
    makepool(file::String, ng::Int, iscircular::Bool = false)

Create a pool of `ng` copies of the genome in the fasta formatted `file`.

!!! note
    The argument `iscircular` is currently not used.
"""
function makepool(file::String, ng::Int = 1, iscircular::Bool = false)
    open(FASTA.Reader, file) do rdr
        makepool(rdr, ng, iscircular)
    end
end

@inline views(u::MoleculePool) = u.seqviews
@inline genome(u::MoleculePool) = u.genome
@inline nmolecules(u::MoleculePool) = length(views(u))

function Base.show(io::IO, p::MoleculePool)
    println(io, "Pool of $(nmolecules(p)) molecules:")
    mx, av, mn = summarize_lengths(views(p))
    if mx == mn
        println(io, " All molecules are of the same size: $mx")
    else
        println(io, " Maximum molecule size: $mx")
        println(io, " Average molecule size: $av")
        println(io, " Minimum molecule size: $mn")
    end
    tags = summarize_tags(views(p))
    if length(tags) > 0
        println(io, " Number of distinct tags: $(length(keys(tags)))")
    end
end

# Transformations
# ---------------

function amplify(p::MoleculePool, n::Int)
    np = MoleculePool(genome(p), amplify(views(p), n))
    return np
end

"""
    fragment(p::MoleculePool, meansize::Int)
    
Create a new pool by breaking up the DNA fragments in an input pool.

This method breaks up a DNA molecule in a pool `p`, such that the average
length of the fragments is approximately `meansize`.

It fragments a molecule by scattering an appropriate number of breakpoints
across the molecule, before cutting the molecule at those breakpoints.

!!! note
    Breakpoints are scattered entirely at random across a molecule.
    No two or more breakpoints can fall in exactly the same place, as those
    positions are sampled without replacement.

!!! note
    The appropriate number of breakpoints to scatter across a molecule is
    calculated as:
    
    ```math
    \\frac{L}{S} - 1
    ```
    
    Where ``L`` is the length of the molecule being fragmented, and ``S`` is
    the desired expected fragment size.
    This calculation assumes breakpoints fall randomly across the molecule
    (see above note).

!!! note
    If a DNA molecule being fragmented is smaller than the desired `meansize`,
    then it will not be broken, it will simply be included in the new pool.
"""
function fragment(p::MoleculePool, meansize::Int)
    np = MoleculePool(genome(p), fragment(views(p), meansize::Int))
    return np
end

"""
    subsample(p::MoleculePool, n::Int)

Create a new pool by sampling an input pool.

!!! note
    DNA molecules in the input pool `p` are selected according to the
    uniform distribution; no one molecule is more or less likely to be selected
    than another.

!!! note
    Sampling is done without replacement, so it is impossible for the new
    pool that is created to recieve one molecule the input pool twice.
"""
function subsample(p::MoleculePool, n::Int)
    np = MoleculePool(genome(p), subsample(views(p), n))
    return np
end

"""
    tag(u::MoleculePool, ntags::Int)

Create a pool of tagged DNA molecules from some input pool.

The new tagged pool has the same DNA molecules as the input pool.
However, each DNA molecule in the new tagged pool will be assigned a tag
in the range of `1:ntags`.

For any tagged molecules in a pool, any other molecules that
are derived from that tagged molecule will inherit the same tag.
For example, if a DNA fragment in a pool is tagged, and then it is
subsequently fragmented during a [`fragment`](@ref) transform, then all the
smaller fragments derived from that long fragment will inherit that long
fragment's tag.

!!! note
    Which fragment gets a certain tag is completely random.
    It is possible for two distinct DNA molecules in a pool to be assigned
    the same tag. The likelihood of that happening depends on the size of the
    tag pool (`ntags`), and the number of fragments in the pool.
"""
function tag(p::MoleculePool, ntags::Int)
    @info "Attaching $(ntags) tags randomly to each molecule in pool..."
    np = MoleculePool(genome(p), tag(views(p), ntags))
    @info "Done"
    return np
end

"""
    flip(p::MoleculePool)

Create a copy of flipped DNA molecules from some input pool.

The flip function lets you randomly flip some of the molecules in a pool to the
opposite orientation they are in.
"""
function flip(p::MoleculePool)
    return MoleculePool(genome(p), flip_views(views(p)))
end

function select(fn::Function, p::MoleculePool)
    np = MoleculePool(genome(p), select(fn, views(p)))
end