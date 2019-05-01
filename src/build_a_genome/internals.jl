
contains(a::UnitRange{Int}, b::UnitRange{Int}) = intersect(a, b) == b

function complimentary(a::UnitRange{Int}, b::UnitRange{Int})
    lower = min(first(a), first(b)):(max(first(a), first(b)) - 1)
    upper = (min(last(a), last(b)) + 1):max(last(a), last(b))
    return lower, upper
end

function pick_available_subregion(region::UnitRange{Int}, size::Int)
    rgnstart = pick_start_point(region, size)
    rgnstop = rgnstart + size - 1
    return rgnstart:rgnstop
end

function pick_start_point(region::UnitRange{Int}, size::Int)
    if length(region) < size
        throw(ArgumentError("region must be longer or equal to size."))
    end
    laststart = last(region) - size + 1
    return rand(first(region):laststart)
end