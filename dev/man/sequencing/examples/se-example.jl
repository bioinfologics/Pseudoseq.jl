using Pseudoseq.Sequencing

m = Molecules("ecoli-ref.fasta")

amp = amplify(5000)

frag = fragment(40000bp)

ssmpl = subsample(50X, 40000bp)

readmaker = makereads()

errmaker = make_substitutions(FixedProbSubstitutions(0.1))

m |> amp |> frag |> ssmpl |> readmaker |> errmaker |> generate("se-reads.fastq")

my_protocol = errmaker ∘ readmaker ∘ ssmpl ∘ frag ∘ amp

m |> my_protocol |> generate("se-reads.fastq")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

