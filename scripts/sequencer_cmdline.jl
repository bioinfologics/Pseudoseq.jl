module SequencerCmdLine

using ArgParse
using SequencingSim
using Distributions
using BioSequences
using ProgressMeter
import BioCore.IO.stream

function parse_cmdline_options()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--sequences"
            help = "The input FASTA file of DNA sequences, assumed to be extracted and amplified to desired quantities."
            arg_type = String
        "--output", "-o"
            help = "Prefix for output FASTQ files."
            arg_type = String
        "--len", "-l"
            help = "Length of reads."
            arg_type = Int
        "--errmean", "-e"
            help = "Mean number of errors per read during sequencing step."
            arg_type = Int
        "--errsd", "-E"
            help = "SD of number of errors per read during sequencing step."
            arg_type = Int
    end
    return parse_args(s, as_symbols = true)
end

function main()
    arguments = parse_cmdline_options()
    errdist = Normal(arguments[:errmean], arguments[:errsd])
    fragments = open(FASTA.Reader, arguments[:sequences])
    fragment_record = FASTA.Record()
    outprefix = arguments[:output]
    r1writer = open(FASTQ.Writer, string(outprefix, "_R1.fastq"))
    r2writer = open(FASTQ.Writer, string(outprefix, "_R2.fastq"))
    
    p = SequencingSim.setup_progbar(fragments, "Fragments sequenced: ")
    
    while !eof(fragments)
        read!(fragments, fragment_record)
        fragment_seq = sequence(SequencingSim.SEQ_TYPE, fragment_record)
        fragname = FASTA.identifier(fragment_record)
        r1read, r2read = SequencingSim.sequence_ends(fragment_seq, arguments[:len], errdist)
        qual = Vector{Int}(undef, length(r1read))
        FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, zeros(length(r1read)), qual, 1, length(r1read))
        r1fastq = FASTQ.Record(string(fragname, "_R1"), r1read, qual)
        FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, zeros(length(r2read)), qual, 1, length(r2read))
        r2fastq = FASTQ.Record(string(fragname, "_R2"), r2read, qual)
        write(r1writer, r1fastq)
        write(r2writer, r2fastq)
        update!(p, position(stream(fragments)))
    end
    close(r1writer)
    close(r2writer)
    close(fragments)
end

end

SequencerCmdLine.main()