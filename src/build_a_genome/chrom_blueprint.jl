struct RepetitionBlueprint
    origin::Int
    destination::Int
    size::Int
end

struct HeterozygosityBlueprint
    position::Int
    alleles::Vector{DNA}
end

position(h::HeterozygosityBlueprint) = h.position
alleles(h::HeterozygosityBlueprint) = h.alleles
Base.isless(a::HeterozygosityBlueprint, b::HeterozygosityBlueprint) = a.position < b.position

struct ChromosomeBlueprint
    seqlen::Int
    available::Set{UnitRange{Int}}
    ncopies::Int
    hets::Vector{HeterozygosityBlueprint}
    intra_repetitions::Vector{RepetitionBlueprint}
end

function ChromosomeBlueprint(len::Int, ncopies::Int)
    return ChromosomeBlueprint(len, Set([1:len]), ncopies, Vector{HeterozygosityBlueprint}(), Vector{RepetitionBlueprint}())
end

function ChromosomeBlueprint(cp::ChromosomeBlueprint, v::Vector{RepetitionBlueprint})
    return ChromosomeBlueprint(cp.seqlen, cp.available, cp.ncopies, cp.hets, v)
end

function ChromosomeBlueprint(cp::ChromosomeBlueprint, v::Vector{HeterozygosityBlueprint})
    return ChromosomeBlueprint(cp.seqlen, cp.available, cp.ncopies, v, cp.intra_repetitions)
end

function ChromosomeBlueprint(cp::ChromosomeBlueprint, s::Set{UnitRange{Int}})
    return ChromosomeBlueprint(cp.seqlen, s, cp.ncopies, cp.hets, cp.intra_repetitions)
end

"Get the number of copies of a chromosome planned in a chromosome blueprint."
ncopies(cp::ChromosomeBlueprint) = cp.ncopies

"Get a copy of the available blocks of sequence in a chromosome blueprint."
available_blocks(cp::ChromosomeBlueprint) = copy(cp.available)

"Get a reference to the available blocks of sequence in a chromosome blueprint."
available_blocks_unsafe(cp::ChromosomeBlueprint) = cp.available

"Get a copy of the repetitions planned in a chromosome blueprint."
repetitions(cp::ChromosomeBlueprint) = copy(cp.intra_repetitions)

"Get a refernece to the repetitions planned in a chromosome blueprint."
repetitions_unsafe(cp::ChromosomeBlueprint) = cp.intra_repetitions

"Get a copy of the heterozygosity planned in a chromosome blueprint."
hets(cp::ChromosomeBlueprint) = copy(cp.hets)

"Get a reference to the heterozygosity planned in a chromosome blueprint."
hets_unsafe(cp::ChromosomeBlueprint) = cp.hets

"Get the length (in bp) of the chromosome planned in a chromosome blueprint."
Base.length(cp::ChromosomeBlueprint) = cp.seqlen

function Base.show(io::IO, cp::ChromosomeBlueprint)
    println(io, "Specification for a chromosome:")
    println(io, " Number of copies: ", ncopies(cp))
    println(io, " Length of chromosomes: ", length(cp))
    ablocks = sort!(collect(available_blocks_unsafe(cp)))
    asum = sum(length(x) for x in ablocks)
    println(io, ' ', asum, "bp available in the following stretches:")
    if length(ablocks) <= 20
        for b in ablocks
            println("  ", first(b), "bp .. ", last(b), "bp")
        end
    else
        for i in 1:10
            b = ablocks[i]
            println("  ", first(b), "bp .. ", last(b), "bp")
        end
        println("    ⋮")
        li = lastindex(ablocks)
        for i in (li - 9):li
            b = ablocks[i]
            println("  ", first(b), "bp .. ", last(b), "bp")
        end
    end
    hlen = length(hets_unsafe(cp))
    if hlen > 0
        println(' ', length(hets_unsafe(cp)), " heterozygous positions at:")
        hs = sort!(hets(cp))
        if hlen <= 20
            for h in hs
                println("  location: ", position(h), "bp  bases: ", string(alleles(h)))
            end
        else
            for i in 1:10
                h = hs[i]
                println("  location: ", position(h), "bp  bases: ", string(alleles(h)))
            end
            println("    ⋮")
            li = lastindex(hs)
            for i in (li - 9):li
                h = hs[i]
                println("  location: ", position(h), "bp  bases: ", string(alleles(h)))
            end
        end
    end
end

function count_het(cb::ChromosomeBlueprint)
    counts = zeros(Int, ncopies(cb), ncopies(cb))
    for het in hets_unsafe(cb)
        for i in 1:ncopies(cb)
            for j in i:ncopies(cb)
                inc = ifelse(het.alleles[i] != het.alleles[j], 1, 0)
                counts[i, j] = counts[j, i] = (counts[i, j] + inc)
            end
        end
    end
    return counts
end

function Base.copy(cp::ChromosomeBlueprint)
    a = available_blocks(cp)
    s = copy(cp.snps)
    r = repetitions(cp)
    return ChromosomeBlueprint(cp.seqlen, a, cp.ncopies, s, r)
end


# Choosing and consuming available regions of the genome

## Internals

function find_big_block(blocks::Set{UnitRange{Int}}, region::UnitRange{Int})
    for block in blocks
        if contains(block, region)
            return block
        end
    end
    throw(ArgumentError(string("Region ", region, " is not contained by an available block.")))
end

function consume_region!(blocks::Set{UnitRange{Int}}, region::UnitRange{Int})
    bigblock = find_big_block(blocks, region)
    pop!(blocks, bigblock)
    l, r = complimentary(bigblock, region)
    if !isempty(l)
        push!(blocks, l)
    end
    if !isempty(r)
        push!(blocks, r)
    end
    return blocks
end


## Public API

"""
    suggest_regions(cp::ChromosomeBlueprint, size::Int, n::Int)
    
A useful utility function to assist planning features in a chromosome blueprint.
    
This method returns a vector of non-overlapping, regions of the chromosome planned
represented by the ChromosomeBlueprint `cb`.

These regions are free regions: they are untouched by any other planned features,
and so may be used when planning other features (See also: [`plan_repetition`](@ref), [`plan_het`](@ref)).

The regions will be `size`bp in length.

!!! warning
    This function was designed for use interactively in a julia session.
    
    If this method cannot find `n` free regions of the size you've asked for, it
    will still give an output vector containing the regions it did manage to find,
    but it will issue a warning to the terminal. This will get increasingly
    likely as you fill the chromosome blueprint up with features.
    
    If you are using this method in a script or program where you depend on a
    reliable number of regions, either add a check to make sure you got the
    number of regions you need, or use this method interactively, and hard code
    an appropriate output into your script.
"""
function suggest_regions(cp::ChromosomeBlueprint, size::Int, n::Int)
    suitable_blocks = filter(x -> length(x) ≥ size, available_blocks_unsafe(cp))
    suggested_regions = Vector{UnitRange{Int}}(undef, n)
    nfound = 0
    while nfound < n
        if isempty(suitable_blocks) # We've run out of big enough available blocks...
            @warn string("Could only suggest ", nfound, " available regions of ", size, "bp or longer")
            break
        end
        bigblock = rand(suitable_blocks)
        pop!(suitable_blocks, bigblock)
        smallblock = pick_available_subregion(bigblock, size)
        l, r = complimentary(bigblock, smallblock)
        if length(l) ≥ size
            push!(suitable_blocks, l)
        end
        if length(r) ≥ size
            push!(suitable_blocks, r)
        end
        suggested_regions[nfound += 1] = smallblock
    end
    return resize!(suggested_regions, nfound)
end



function consume_regions(cp::ChromosomeBlueprint, blocks::Vararg{UnitRange{Int}})
    avb = available_blocks(cp)
    for b in blocks
        consume_region!(avb, b)
    end
    return ChromosomeBlueprint(cp, avb)
end

"Returns true if the `region` of the chromosome in the blueprint is unoccupied by a planned feature."
function isavailable(cp::ChromosomeBlueprint, region::UnitRange{Int})
    for block in available_blocks_unsafe(cp)
        contains(block, region) && return true
    end
    return false
end

function checkavailability(cp::ChromosomeBlueprint, from::UnitRange{Int}, to::UnitRange{Int})
    if !isavailable(cp, from)
        throw(ArgumentError("The first region is not available"))
    end
    if !isavailable(cp, to)
        throw(ArgumentError("The second region is not available"))
    end
end

function checkavailability(cp::ChromosomeBlueprint, interval::UnitRange{Int})
    if !isavailable(cp, interval)
        throw(ArgumentError("requested region not available: $interval"))
    end
end







