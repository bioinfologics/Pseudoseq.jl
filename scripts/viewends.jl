module Viewends

using ArgParse
using SequencingSim

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input", "-I"
            help = "The path of an input workspace file"
            arg_type = String
        "--output", "-O"
            help = "Prefix for output workspace file"
            arg_type = String
        "--len", "-L"
            help = "Length of the ends"
            arg_type = Int
    end
    settings = parse_args(ARGS, s, as_symbols = true)
    vs = load_views(settings[:input]::String)
    newvs = take_paired_ends(vs, settings[:len]::Int)
    save_views(newvs, settings[:output]::String)
    return 0
end

end # module

Viewends.julia_main(ARGS)