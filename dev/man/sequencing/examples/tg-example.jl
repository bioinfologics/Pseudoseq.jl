using Pseudoseq.Sequencing

m = Molecules("ecoli-ref.fasta")

amp = amplify(5000)

firstfrag = fragment(40000bp)

tagger = tag(1000000)

secondfrag = fragment(700bp)

size_filter = select(x -> 900 >= length(x) >= 450)

ssmpl = subsample(50X, 2 * 250bp)

readmaker = makereads(2 * 250bp)

errmaker = make_substitutions(FixedProbSubstitutions(0.001))

m |> amp |> firstfrag |> tagger |> secondfrag |> size_filter |> ssmpl |> readmaker |> errmaker |> generate("tagged_reads.fastq")

my_protocol = errmaker ∘ readmaker ∘ ssmpl ∘ size_filter ∘ secondfrag ∘ tagger ∘ firstfrag ∘ amp

m |> my_protocol |> generate("tagged_reads.fastq")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

