using Pseudoseq, Pseudoseq.Sequencing, Documenter, Literate

Documenter.post_status(; type = "pending", repo = "github.com/bioinfologics/Pseudoseq.jl.git")

SEQ_EXAMPLES = [joinpath("examples/sequencing", f) for f in ("pe-example.jl",
                                                             "se-example.jl",
                                                             "tg-example.jl")]

#BAG_EXAMPLES = [joinpath("examples/build-a-genome", f) for f in ("build-a-yeast.jl",)]


OUTPUT = joinpath(@__DIR__, "src/man")
SEQ_OUTPUT = joinpath(OUTPUT, "sequencing/examples")
#BAG_OUTPUT = joinpath(OUTPUT, "build-a-genome/examples")
mkdir(SEQ_OUTPUT)
#mkdir(BAG_OUTPUT)

cp("examples/sequencing/ecoli-ref.fasta", "docs/src/man/sequencing/examples/ecoli-ref.fasta")
#cp("examples/build-a-genome/yeast-chr1.fasta", "docs/src/man/build-a-genome/examples/yeast-chr1.fasta")

for ex in SEQ_EXAMPLES
    Literate.markdown(ex, SEQ_OUTPUT)
    Literate.notebook(ex, SEQ_OUTPUT)
    Literate.script(ex, SEQ_OUTPUT)
end

#for ex in BAG_EXAMPLES
#    Literate.markdown(ex, BAG_OUTPUT)
#    Literate.notebook(ex, BAG_OUTPUT)
#    Literate.script(ex, BAG_OUTPUT)
#end

makedocs(
    format = Documenter.HTML(
        prettyurls = haskey(ENV, "GITHUB_ACTIONS")
    ),
    modules = [Pseudoseq, Sequencing],
    sitename = "Pseudoseq.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "PuzzleMaker" => "man/puzzlemaker/concepts.md",
            "Sequencing" => [
                "Core concepts & workflow" => "man/sequencing/concepts.md",
                "Examples" => [
                    "Paired end reads" => "man/sequencing/examples/pe-example.md",
                    "Long single end reads" => "man/sequencing/examples/se-example.md",
                    "Tagged paired end reads" => "man/sequencing/examples/tg-example.md"
                ]
            ]
        ],
        "API" => [
            "Sequencing" => [
                "Molecule Pool" => "api/pool.md",
                "Reads" => "api/reads.md",
            ]
        ]
    ],
    authors = "Ben J. Ward."
)

# Scripts will generate some read file output FASTQs, we don't really want to
# keep these, and their size might trip up CI when we upload docs to gh-pages,
# so we find all files ending in .fastq, and delete them.
run(`find . -name "*.fastq" -type f -delete`)

deploydocs(
    repo = "github.com/bioinfologics/Pseudoseq.jl.git",
    push_preview = true
)