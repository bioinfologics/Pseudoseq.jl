using Pseudoseq.Sequencing

m = Molecules("ecoli-ref.fasta")

amp = amplify(5000)

frag = fragment(700bp)

size_filter = select(x -> 900 >= length(x) >= 450)

ssmpl = subsample(50X, 2 * 250bp)

readmaker = makereads(2 * 250bp)

errmaker = make_substitutions(FixedProbSubstitutions(0.001))

m |> amp |> frag |> size_filter |> ssmpl |> readmaker |> errmaker |> generate("pe-reads.fastq")

my_protocol = errmaker ∘ readmaker ∘ ssmpl ∘ size_filter ∘ frag ∘ amp

m |> my_protocol |> generate("pe-reads.fastq")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

