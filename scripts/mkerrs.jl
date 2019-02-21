module Mkerrs

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
        "--err", "-R"
            help = "The per base error rate"
            arg_type = Float64
    end
    settings = parse_args(ARGS, s, as_symbols = true)
    vs = load_views(settings[:input]::String)
    errs = make_errors(vs, settings[:err]::Float64)
    save_views_and_errs(vs, errs, settings[:output]::String)
    return 0
end

end # module

Mkerrs.julia_main(ARGS)