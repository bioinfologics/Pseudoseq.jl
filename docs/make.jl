using Documenter, Literate, Pseudoseq

SEQ_EXAMPLES = [joinpath("examples/sequencing", f) for f in ("pe-example.jl",
                                                             "se-example.jl",
                                                             "tg-example.jl")]

BAG_EXAMPLES = [joinpath("examples/build-a-genome", f) for f in ("build-a-yeast.jl",)]


OUTPUT = joinpath(@__DIR__, "src/examples")
SEQ_OUTPUT = joinpath(OUTPUT, "sequencing")
BAG_OUTPUT = joinpath(OUTPUT, "build-a-genome")
mkdir(OUTPUT)
mkdir(SEQ_OUTPUT)
mkdir(BAG_OUTPUT)

cp("examples/sequencing/ecoli-ref.fasta", "docs/src/examples/sequencing/ecoli-ref.fasta")
cp("examples/build-a-genome/yeast-chr1.fasta", "docs/src/examples/build-a-genome/yeast-chr1.fasta")

for ex in SEQ_EXAMPLES
    Literate.markdown(ex, SEQ_OUTPUT)
    Literate.notebook(ex, SEQ_OUTPUT)
    Literate.script(ex, SEQ_OUTPUT)
end

for ex in BAG_EXAMPLES
    Literate.markdown(ex, BAG_OUTPUT)
    Literate.notebook(ex, BAG_OUTPUT)
    Literate.script(ex, BAG_OUTPUT)
end

makedocs(
    format = :html,
    modules = [Pseudoseq],
    sitename = "Pseudoseq.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Sequencing" => [
            "Core concepts & workflow" => "sequencing.md",
            "Examples" => [
                "Paired end reads" => "examples/sequencing/pe-example.md",
                "Long single end reads" => "examples/sequencing/se-example.md",
                "Tagged paired end reads" => "examples/sequencing/tg-example.md"
            ]
        ],
        "Build-a-Genome" => [
            "Core concepts & workflow" => "build-a-genome.md",
            "Examples" => [
                "Build-a-Yeast" => "examples/build-a-genome/build-a-yeast.md",
            ]
        ],
        
        "API" => [
            "Sequencing" => [
                "sequence" => "api/sequence.md",
                "Molecule Pool" => "api/pool.md",
                "Reads" => "api/reads.md",
            ]
            "Build-a-Genome" => "api/chromosome-blueprint.md"
        ]
    ],
    authors = "Ben J. Ward."
)

# Scripts will generate some read file output FASTQs, we don't really want to
# keep these, and their size might trip up CI when we upload docs to gh-pages,
# so we find all files ending in .fastq, and delete them.
run(`find . -name "*.fastq" -type f -delete`)

deploydocs(
    repo = "github.com/bioinfologics/Pseudoseq.jl.git"
)