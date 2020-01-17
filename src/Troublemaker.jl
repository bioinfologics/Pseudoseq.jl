module Troublemaker

export
    MotifStitcher,
    RandomMotif,
    add_motif!,
    add_motifs!,
    add_motif_arrangement!,
    generate_random_motif_sequences!,
    stitch_motifs


using BioSequences
using GenomeGraphs

const DEFAULT_SAMPLER = SamplerWeighted(dna"ATCG", fill(0.25, 3))
const SEQ_T = LongSequence{DNAAlphabet{2}}

struct RandomMotif
    len::Int
    sampler::SamplerWeighted{DNA}
end

struct SiblingMotif
    source_motif::Int
    poly_positions::Int
end

mutable struct MotifStitcher
    n_motifs::Int
    random_motifs::Dict{Int,RandomMotif}
    fixed_motifs::Dict{Int,SEQ_T}
    motif_similarities::Dict{Int,SiblingMotif}
    motif_order::Vector{Vector{Int}}
end

function MotifStitcher()
    return MotifStitcher(
        0,
        Dict{Int,RandomMotif}(),
        Dict{Int,SEQ_T}(),
        Dict{Int,SiblingMotif}(),
        Vector{Vector{Int}}(),
    )
end

@inline function add_motif!(ms::MotifStitcher, rm::RandomMotif)
    ms.n_motifs += 1
    ms.random_motifs[ms.n_motifs] = rm
    return ms
end

@inline function add_motif!(ms::MotifStitcher, len::Int)
    ms.n_motifs += 1
    ms.random_motifs[ms.n_motifs] = RandomMotif(len, DEFAULT_SAMPLER)
    return ms
end

@inline function add_motif!(ms::MotifStitcher, fm::SEQ_T)
    ms.n_motifs += 1
    ms.fixed_motifs[ms.n_motifs] = fm
    return ms
end

function add_motifs!(ms::MotifStitcher, args...)
    for arg in args
        add_motif!(ms, arg)
    end
    return ms
end

function add_motif_arrangement!(ms::MotifStitcher, arr::Vector{Int})
    arr′ = abs.(arr)
    @assert minimum(arr′) === 1
    @assert maximum(arr′) ≤ ms.n_motifs
    push!(ms.motif_order, copy(arr))
    return ms
end









#=
function fix_motif_sequence!(ms::MotifStitcher, motifid::Int, seq::SEQ_T)
    delete!(ms.random_motifs, motifid)
    ms.fixed_motifs[motifid] = seq
    ms.motif_lengths[motifid] = length(seq)
    return ms
end
=#

#=
function _delete_and_rekey_dictionary!(d::Dict{Int,T}, key::Int) where {T}
    rekeyed = Dict(k-1 => v for (k, v) in d if k > key)
    for k in keys(d)
        if k >= key
            delete!(d, k)
        end
    end
    merge!(d, rekeyed)
    return d
end

function delete_motif!(ms::MotifStitcher, motifid::Int)
    deleteat!(ms.motif_lens, motifid)
    _delete_and_rekey_dictionary!(ms.random_motifs, motifid)
    _delete_and_rekey_dictionary!(ms.fixed_motifs, motifid)
    _delete_and_rekey_dictionary!(ms.motif_similarities, motifid)
    # TODO: Check the motif similarities, if any of them are based off the deleted motif, 
    # issue a warning and convert it to a random sequence.
    return ms
end
=#

#=
function set_motif_lengths!(pm::MotifStitcher, lens::Vector{Int})
    @assert all(lens .> 0)
    resize!(pm.motif_lens, length(lens))
    copyto!(pm.motif_lens, lens)
    return pm
end
=#



function generate_motif_sequences(ms::MotifStitcher)
    sequences = Dict{Int, SEQ_T}()
    for (k, v) in ms.random_motifs
        sequences[k] = randseq(DNAAlphabet{2}(), v.sampler, v.len)
    end
    merge!(sequences, ms.fixed_motifs)
    # TODO: Generate sibling motifs!
    return sequences
end

function stitch_motifs(ms::MotifStitcher)
    sequences = Vector{SEQ_T}(undef, length(ms.motif_order))
    motif_seqs = generate_motif_sequences(ms)
    for i in eachindex(ms.motif_order)
        s = SEQ_T()
        mo = ms.motif_order[i]
        for m in mo
            if m < 0
                append!(s, reverse_complement(motif_seqs[abs(m)]))
            else
                append!(s, motif_seqs[m])
            end
        end
        sequences[i] = s
    end
    return sequences
end

function GenomeGraphs.Graphs.SequenceDistanceGraph(ms::MotifStitcher)
    sequences = generate_motif_sequences(ms)
    sdg = GenomeGraphs.Graphs.SequenceDistanceGraph{SEQ_T}()
    for i in 1:length(sequences)
        GenomeGraphs.Graphs.add_node!(sdg, sequences[i])
    end
    for path in ms.motif_order
        for i in 2:lastindex(path)
            GenomeGraphs.Graphs.add_link!(sdg, -path[i - 1], path[i], 0)
        end
    end
    return sdg
end

end # module