using Documenter, Pseudoseq

makedocs(
    format = :html,
    modules = [Pseudoseq],
    sitename = "Pseudoseq.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "User manual" => [
            "Creating the universe" => "man/create_universe.md",
            "Fragment molecules" => "man/fragment_molecules.md",
            "Subsample molecules" => "man/subsample_molecules.md",
            "Generate reads" => "man/generate_reads.md"
        ],
        "Developer info" =>  [
            "Views" => "dev/views.md",
        ]
    ],
    authors = "Ben J. Ward."
)

deploydocs(
    repo = "github.com/bioinfologics/Pseudoseq.jl.git",
    deps = nothing,
    make = nothing
)