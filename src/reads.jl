
abstract type Reads end

abstract type PairedReads <: Reads end

struct PairedEndReads <: PairedReads
    seqviews::Views
end

struct TenXReads <: PairedReads
    seqviews::Views
end

struct LongReads <: Reads
    seqviews::Views
end



# Paired end specific
function take_ends(sv::SequencingView, length::Int)
    return subview(sv, firstindex(sv), length), subview(sv, lastindex(sv), lastindex(sv) - length + 1)
end

function make_paired_reads(vs::Views, len::Int)
    i = 0
    newviews = Views(length(vs) * 2)
    for view in vs
        if length(view) > len # Guard against sequences that are too small.
            fwdview, revview = take_paired_ends(view, len)
            newviews[i += 1] = fwdview
            newviews[i += 1] = revview
        end
    end
    resize!(newviews, i)
    return newviews
end

function make_errors(vs::Views, rate::Float64)
    nbases = sum([length(v) for v in vs])
    nerrs = Int(floor(nbases * rate))
    errs = Dict{Int, Vector{Int}}()
    reads = rand(1:length(vs), nerrs)
    positions = [rand(1:length(vs[read])) for read in reads]
    return reads, positions
end

# User facing method.
function make_reads(vs::Views, len::Int, err::Float64, mode::Symbol)
    if mode === :pe
        return make_paired_reads(vs, len)
    elseif mode === :linked
        
    elseif mode === :long
        
    end
end




function extract_sequence(ref::Sequence, view::SequencingView)
    if isbwd(view)
        subseq = ref[last(view):first(view)]
        reverse_complement!(subseq)
    else
        subseq = ref[first(view):last(view)]
    end
    return subseq
end

const nucs = [ACGT...]

function add_errors!(seq::Sequence, errs::Vector{Int})
    for err in errs
        truebase = seq[err]
        posbases = filter(!isequal(truebase), nucs)
        seq[err] = rand(posbases)
    end
end