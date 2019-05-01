# Repeated content
# ----------------

"""
    plan_repetition(cb::ChromosomeBlueprint, from::UnitRange{Int}, to::UnitRange{Int})
    
Plan a repetition in a chromosome, where the bases in the `from` region of the
chromosome, are guaranteed to occur again in the `to` region of the chromosome.

Every copy of the chromosome will have this same repetition.

Creates a new chromosome blueprint, based on the input blueprint `cb`.

!!! note
    In the new blueprint, the `to` region of the planned chromosome
    will be consumed, and cannot be used to plan any other subsequently added
    features.
"""
function plan_repetition(cb::ChromosomeBlueprint, from::UnitRange{Int}, to::UnitRange{Int})
    if length(from) != length(to)
        throw(ArgumentError("from and to must be the same length"))
    end
    checkavailability(cb, to)
    newspec = consume_regions(cb, to)
    newrep = RepetitionBlueprint(first(from), first(to), length(from))
    return ChromosomeBlueprint(newspec, push!(repetitions(newspec), newrep))
end


"""
    plan_repetition(cb::ChromosomeBlueprint, from::Int, to::Int, size::Int)
    
Plan a repetition in a chromosome, where the bases in the
`from`:(`from` + `size` - 1) region of the chromosome, are guaranteed to occur
again in the `to`:(`to` + `size` - 1) region of the chromosome.

Every copy of the chromosome will have this same repetition.

Creates a new chromosome blueprint, based on the input blueprint `cb`.
"""
function plan_repetition(cb::ChromosomeBlueprint, from::Int, to::Int, size::Int)
    rfrom = from:(from + size - 1)
    rto = to:(to + size - 1)
    return plan_repetition(cb, rfrom, rto)
end


"""
    plan_repetition(cb::ChromosomeBlueprint, intervals::Vector{UnitRange{Int}})
    
A conveinience method of `plan_repetition`.
Designed to ease the process of planing a series of repetitions in a chromosome.

Every copy of the chromosome will have these same repetitions.

Creates a new chromosome blueprint, based on the input blueprint `cb`.

!!! tip
    Use the [`suggest_regions`](@ref) function to help decide on a set of
    sites to make heterozygous.

!!! note
    The number of intervals provided must be an even number. This is because
    intervals 1 & 2  define the first repeat, intervals 3 & 4 define the second, and
    so on.

!!! note
    In the new blueprint, the regions the repetitions occupy, have been consumed,
    and cannot be used to plan any other subsequently added features.
"""
function plan_repetition(cb::ChromosomeBlueprint, intervals::Vector{UnitRange{Int}})
    @assert iseven(length(intervals))
    lastcb = cb
    for i in firstindex(intervals):2:lastindex(intervals)
        lastcb = plan_repetition(lastcb, intervals[i], intervals[i + 1])
    end
    return lastcb
end