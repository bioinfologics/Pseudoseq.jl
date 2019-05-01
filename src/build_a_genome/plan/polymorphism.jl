# Heterozygous content
# --------------------

### The core heterozygosity planning method. All others eventually call this one.

"""
    plan_het(cb::ChromosomeBlueprint, pos::Int, alleles::Vector{DNA})

Plan heterozygosity between copies of a chromosome at position `pos`.

The `alleles` vector must contain a nucleotide for each copy of the
chromosome.

For example if you provided a vector of `[DNA_A, DNA_C]`, for the heterozygous
site you're defining at `pos`, the first copy of the chromosome will have an `A`,
and the second copy of the chromosome will have a `C`.

Creates a new chromosome blueprint, based on the input blueprint `cb`.

!!! tip
    Use the [`suggest_alleles`](@ref) function to help decide on a set of
    alleles to use.

!!! note
    In the new blueprint, the positions that were used to plan the heterozygosity,
    have been consumed, and cannot be used to plan any other subsequently added
    features.
"""
function plan_het(cb::ChromosomeBlueprint, pos::Int, alleles::Vector{DNA})
    if length(alleles) != ncopies(cb)
        throw(ArgumentError(string("alleles must have ", ncopies(cb), " entries")))
    end
    for a in alleles
        if isambiguous(a) 
            throw(ArgumentError("alleles must not contain any ambiguous nucleotides"))
        end
    end
    if length(unique(alleles)) == 1
        throw(ArgumentError("alleles must contain at least two different nucleotides"))
    end
    checkavailability(cb, pos:pos)
    interplan = consume_regions(cb, pos:pos)
    newhet = HeterozygosityBlueprint(pos, alleles)
    return ChromosomeBlueprint(interplan, push!(hets(interplan), newhet))
end

### Allele suggestion methods.

"""
    suggest_alleles(groups::Vector{Int})

Suggest an allele pattern for a single heterozygous site.

The `groups` vector should have a number of 
The nth value in the `groups` designates a group to the nth chromosome copy.

E.g. consider for a triploid, the vector `[2, 1, 2]` means that there are two
groups (group 1 and group 2), and that the second chromosome copy of the three
belongs to group 1, and the other two copies belong to group 2.

Chromosome copies that belong to the same group will get the same base at each
heterozygous site.

For the vector `[2, 1, 2]`, this method might return a suggested allele pattern
such as `[DNA_A, DNA_G, DNA_A]`.

!!! note
    Which base corresponds to which group number will be determined at random.
"""
function suggest_alleles(groups::Vector{Int})
    @assert all(g ≤ 4 for g in groups)
    possible = Random.shuffle!([DNA_A, DNA_C, DNA_G, DNA_T])
    return possible[groups]
end

"""
    suggest_alleles(groups::Vector{Int}...)

Suggest an allele pattern for a single heterozygous site.

Multiple integer vectors make up `groups`.

Each vector defines a group, and each value in a vector determines which copies
are allocated to that group.

E.g. consider for a triploid, two vectors of `[2] and [1, 3]` passed as arguments
to `groups`. They define two groups, the first vector contains a number
denoting the second chromosome copy, and he second vector contains numbers
representing the first and third chromosome copy.

For the vectors `[2]` and `[1, 3]`, this method might return a suggested allele
pattern such as `[DNA_A, DNA_G, DNA_A]`.

!!! note
    Which base corresponds to which group number will be determined at random.
"""
function suggest_alleles(groups::Vector{Int}...)
    l = sum(length(g) for g in groups)
    newgroups = Vector{Int}(undef, l)
    for (i, g) in enumerate(groups)
        for j in g
            newgroups[j] = i
        end
    end
    return suggest_alleles(newgroups)
end

"""
    suggest_alleles(ncopies::Int, ngroups::Int)

Suggest an allele pattern for a single heterozygous site.

By providing a number representing the number of copies of a chromosome,
and a number representing the number of groups (which at a minimum must be 2).
This method of `suggest_alleles` will randomly allocate the `ncopies` of the
chromosome into the `ngroups` groups.

For example, if you ask for a suggested allele pattern for 3 copies of a chromosome,
and 2 groups, this method might return a suggested pattern of `[DNA_G, DNA_A, DNA_G]`.

!!! note
    Which base corresponds to which group number will be determined at random.
"""
function suggest_alleles(ncopies::Int, ngroups::Int)
    @assert ngroups ≤ ncopies
    @assert ngroups > 1
    groups = Random.shuffle!(collect(Iterators.take(Iterators.cycle(Random.shuffle!(collect(Base.OneTo(ngroups)))), ncopies)))
    return suggest_alleles(groups)
end

"""
    suggest_alleles(npositions::Int, args...)

Suggest allele patterns for numerous heterozygous sites.

This method of `suggest_alleles` calls another `suggest_alleles` with `args...`
`npositions` times, and returns a vector of the results of each call.
"""
function suggest_alleles(npositions::Int, args...)
    return [suggest_alleles(args...) for _ in 1:npositions]
end

### Position and allele argument preparations, used to adapt the generic `plan_het` method's behaviour.

_prepare_het_positions(cb::ChromosomeBlueprint, pos::Int) = _prepare_het_positions(cb, suggest_regions(cb, 1, pos))
_prepare_het_positions(cb::ChromosomeBlueprint, pos::Float64) = _prepare_het_positions(cb, Int(floor(length(cb) * pos)))
_prepare_het_positions(cb::ChromosomeBlueprint, pos::Vector{Int}) = pos 
function _prepare_het_positions(cb::ChromosomeBlueprint, pos::Vector{UnitRange{Int}})
    for p in pos
        @assert length(p) == 1
    end
    return [first(p) for p in pos]
end

_prepare_alleles(cb::ChromosomeBlueprint, npos::Int, x::Vector{Vector{DNA}}) = x
_prepare_alleles(cb::ChromosomeBlueprint, npos::Int, ngroups::Int) = suggest_alleles(npos, ncopies(cb), ngroups)
_prepare_alleles(cb::ChromosomeBlueprint, npos::Int, args...) = suggest_alleles(npos, args...)

"""
    plan_het(cb::ChromosomeBlueprint, pos, alle...)

A generic method of `plan_het`.

Creates a new chromosome blueprint, based on the input blueprint `cb`.

The `pos` argument determines which sites in a chromosome are made heterozygous,
and depending on the type of argument provided, the behaviour differs slightly:

- If `pos` is an integer, then that many available sites are selected at random
  to be heterozygous.

- If `pos` is a float, then it is treated as a proportion of the length of the,
  chromosome, which is converted into an integer, and that many sites are
  selected at random to be heterozygous.

- If `pos` is a vector of integers, it is trated as a list of positions the user
  wants to be heterozygous.

The `alle...` argument is a vararg argument which defines the how bases will be 
allocated at each heterozygous site.

Depending on the type of the argument provided, the behaviour differs slightly,
based on the different methods of [`suggest_alleles`](@ref):

- If `alle...` is a single integer value. It is taken to mean the number of states
  each heterozygous site has. At a minumum this number must be two, as a site with
  only one state is not heterozygous by definition. For each site, which chromosome
  copies get which state is determined at random.

- If `alle...` is a vector of integers (one integer for every chromosome copy),
  then the nth value in the vector dictates the group the nth chromosome copy
  belongs to. For example for a triploid, the vector `[2, 1, 2]` means that
  there are two groups (group 1 and group 2), and that the second chromosome
  copy of the three belongs to group 1, and the other two copies belong to group
  2. Chromosome copies that belong to the same group will get the same base at
  each heterozygous site. In our example, copies 1 and 3 will get the same base,
  as they have been given the same group number.
  Chromosome copy 2 will end up with a different base. Which base corresponds to
  which group number will be determined at random.

- If `alle...` is multiple vectors of integers, then each vector constitutes a
  group, and each value in the vector determines which copies are in that group.
  It is alternative way to giving the same information as you can with a
  single vector (as above). E.g. consider the example of `[2, 1, 2]` from the
  previous point. This can be expressed as two vectors of `[2] and [1, 3]`:
  The first vector contains a number indicating the second chromosome copy, and
  the second vector contains numbers representing the first and third chromosome
  copy.

- If `alle...` is a vector of nucleotide vectors, then the nth vector of
  nucleotides determines which chromosome copy recieves which base at the nth
  heterozygous position. For example, for a triploid, `[[DNA_A, DNA_A, DNA_T], [DNA_G, DNA_C, DNA_C]]`
  would mean that at the first heterozygous site, the first and second copies
  of the chromosome would recieve an `DNA_A` base, and the third copy would
  recieve a `DNA_T` base. For the second heterozygous position, the first copy
  of the chromosome would recieve a DNA_G base, and the other two copies will
  recieve a `DNA_C` base.

(See also: [`plan_repetition`](@ref))
"""
function plan_het(cb::ChromosomeBlueprint, pos, alle...)
    positions = _prepare_het_positions(cb, pos)
    for p in positions
        if !(1 ≤ p ≤ length(cb))
            throw(ArgumentError(string("Position ", p, " is outside of the limits of the chromosome (", 1:length(cb), ')')))
        end
    end
    alleles = _prepare_alleles(cb, length(positions), alle...)
    
    @assert length(positions) == length(alleles)
    lastcb = cb
    for (p, a) in zip(positions, alleles)
        lastcb = plan_het(lastcb, p, a)
    end
    return lastcb
end