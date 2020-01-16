module Troublemaker

export
    MotifStitcher,
    set_motif_lengths!,
    n_motifs,
    set_motif_arrangement!,
    generate_random_motif_sequences!,
    

using BioSequences
using GenomeGraphs

struct MotifStitcher
    motif_lens::Vector{Int}
    motif_order::Vector{Int}
    motif_sequences::Vector{LongSequence{DNAAlphabet{2}}}
end

MotifStitcher() = MotifStitcher(Int[], Int[], LongSequence{DNAAlphabet{2}}[])

function MotifStitcher(lens::Vector{Int}, arr::Vector{Int})
    pm = set_motif_arrangement!(set_motif_lengths!(MotifStitcher(), lens), arr)
end

function set_motif_lengths!(pm::MotifStitcher, lens::Vector{Int})
    @assert all(lens .> 0)
    resize!(pm.motif_lens, length(lens))
    copyto!(pm.motif_lens, lens)
    return pm
end

@inline n_motifs(pm::MotifStitcher) = length(pm.motif_lens)

function set_motif_arrangement!(pm::MotifStitcher, arr::Vector{Int})
    arr′ = abs.(arr)
    @assert minimum(arr′) === 1
    @assert maximum(arr′) ≤ n_motifs(pm)
    resize!(pm.motif_order, length(arr))
    copyto!(pm.motif_order, arr)
    return pm
end

function generate_random_motif_sequences!(pm::MotifStitcher)
    resize!(pm.motif_sequences, length(pm.motif_lens))
    @inbounds for (i, len) in enumerate(pm.motif_lens)
        pm.motif_sequences[i] = randseq(DNAAlphabet{2}(), len)
    end
    return pm
end

function stitch_motifs(pm::MotifStitcher)
    outseq = LongSequence{DNAAlphabet{2}}()
    @inbounds for x in pm.motif_order
        if x < 0
            append!(outseq, reverse_complement(pm.motif_sequences[abs(x)]))
        else
            append!(outseq, pm.motif_sequences[x])
        end
    end
    return outseq
end

function generate(pm::MotifStitcher)
    generate_random_motif_sequences!(pm)
    return stitch_motifs(pm)
end

function GenomeGraphs.Graphs.SequenceDistanceGraph(ms::MotifStitcher)
    sdg = GenomeGraphs.Graphs.SequenceDistanceGraph{eltype(ms.motif_sequences)}()
    for motif in ms.motif_sequences
        GenomeGraphs.Graphs.add_node!(sdg, motif)
    end
    path = ms.motif_order
    for i in 2:lastindex(path)
        GenomeGraphs.Graphs.add_link!(sdg, -path[i - 1], path[i], 0)
    end
    return sdg
end

end # module