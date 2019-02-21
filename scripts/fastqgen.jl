module FastqGen

using ArgParse
using SequencingSim
using BioSequences

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    s = ArgParseSettings()
    @add_arg_table s begin
        "--reference", "-R"
            help = "The path of an input reference genome"
            arg_type = String
            required = true
        "--workspace", "-W"
            help = "The path of an input workspace file"
            arg_type = String
            required = true
        "--output", "-O"
            help = "Prefix for output FASTQ file"
            arg_type = String
            required = true
        "--mode", "-M"
            help = "The mode the generator should use"
            arg_type = String
            default = "pe"
            range_tester = x -> x ∈ ["pe", "10x", "long"]
    end
    settings = parse_args(ARGS, s, as_symbols = true)
    mode = settings[:mode]::String
    
    # Load the views and their errors...
    vs, errs = load_views_and_errors(settings[:workspace]::String)
    refseqs = open(FASTA.Reader, settings[:reference]::String) do rdr
        sequence.(BioSequence{DNAAlphabet{2}}, rdr)
    end
    open(FASTQ.Writer, settings[:output]::String) do wtr
        for i in eachindex(vs)
            v = vs[i]
            e = get(errs, i, Int[])
            ref = refseqs[seqid(v)]
            
            # Extract subsequence from the reference and add errors.
            seq = extract_sequence(ref, v)
            add_errors!(seq, e)
            
            if mode == "10x"
                mer = BioSequence{DNAAlphabet{2}}(DNAKmer{16}(tag(v)))
                seq = mer * seq
            end
            
            # Make the quality string for the read, (currently no real meaning!).
            qual = Vector{Int}(undef, length(seq))
            FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, fill(30, length(seq)), qual, 1, length(seq))
            
            # Create the name for the read...
            fragname = string("Refseq_", seqid(v))
            if mode != "long"
                fragname = string(fragname, ifelse(isodd(i), "_R1", "_R2"))
            end
            
            fqread = FASTQ.Record(fragname, seq, qual)
            write(wtr, fqread)
        end
    end
    
    return 0
end

end # module

FastqGen.julia_main(ARGS)