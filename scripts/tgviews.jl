module Tgviews

using ArgParse
using SequencingSim

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input", "-I"
            help = "The input workspace file"
            arg_type = String
        "--output", "-O"
            help = "Prefix for output workspace file"
            arg_type = String
        "--numtags", "-N"
            help = "The number of tags that exist"
            arg_type = Int
    end
    settings = parse_args(ARGS, s, as_symbols = true)
    views = load_views(settings[:input]::String)
    newviews = tag_views(views, settings[:numtags]::Int)
    save_views(newviews, settings[:output]::String)
    return 0
end

end # module

Tgviews.julia_main(ARGS)