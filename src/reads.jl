# Paired end reads
# ================

abstract type Technology end
abstract type PairedReads <: Technology end
abstract type SingleReads <: Technology end

struct PairedEnd <: PairedReads
    flen::Int
    rlen::Int
end

struct SingleEnd{L<:Union{UInt,Nothing}} <: SingleReads
    len::L
end

struct TaggedPairs <: PairedReads
    flen::Int
    rlen::Int
end

struct Reads{T <: Technology}
    tech::T
    genome::Vector{LongSequence{DNAAlphabet{2}}}
    seqviews::Views
    errpos::Dict{Int, Vector{Int}}
end

@inline views(x::Reads) = x.seqviews
@inline genome(x::Reads) = x.genome
@inline errors(x::Reads) = x.errpos
@inline nreads(x::Reads) = length(x.seqviews)
@inline pick_read(reads::Reads) = rand(1:nreads(reads))
@inline nerrors(x::Reads) = length(x.errpos) > 0 ? sum(length(x) for x in values(x.errpos)) : 0

function total_bp(x::Reads{<:PairedReads})
    halfn = div(nreads(x), 2) 
    return (halfn * x.tech.flen) + (halfn * x.tech.rlen)
end

function Base.show(io::IO, reads::Reads)
    println(summary(reads))
    println(printlen(reads))
    println(string(" Number of errors: ", nerrors(reads)))
end

function mark_errors!(reads::Reads, rate::Float64)
    nbases = total_bp(reads)
    nerrs = Int(floor(nbases * rate))
    pick_error_positions!(reads, nerrs)
end

"""
    mark_errors(reads::Reads, rate::Float64)
    
Create a new set of reads, with errors, from an input set of reads.

When you first create a set of reads using a [`make_reads`](@ref) method,
all the reads in that set are perfect reads. In real sequencing experiments,
sequencers make errors when reading a DNA molecule, and are characterised by
an error rate. This function lets you simulate this characteristic of
sequencers, by marking (at random) positions in the reads that are destined to
be errors in the output FASTQ.

!!! note
    Currently, every position in every read is equally likely to be marked as an
    error.

!!! note
    The number of errors introduced into the reads can be calculated.
    ```math
    E = NR
    ```
    Where ``E`` is the number of errors, ``N`` is the total number of bases
    in your set of reads, and ``R`` is the error rate you provide to this
    function.
"""
function mark_errors(reads::Reads, rate::Float64)
    newreads = typeof(reads)(reads.tech, reads.genome, reads.seqviews, Dict{Int, Vector{Int}}())
    mark_errors!(newreads, rate)
    return newreads
end

function pick_error_positions!(reads::Reads{<:PairedReads}, nerrs::Int)
    errs = errors(reads)
    empty!(errs)
    flength = reads.tech.flen
    rlength = reads.tech.rlen
    for i in 1:nerrs
        readidx = rand(1:nreads(reads))
        pos = rand(1:ifelse(iseven(readidx), rlength, flength))
        if haskey(errs, readidx)
            push!(errs[readidx], pos)
        else
            errs[readidx] = [pos]
        end
    end
end

const nucs = [ACGT...]

function add_errors!(seq::Sequence, errs::Vector{Int})
    for err in errs
        truebase = seq[err]
        posbases = filter(!isequal(truebase), nucs)
        seq[err] = rand(posbases)
    end
end

"""
    generate(filename::String, reads::Reads)
    
Write the `reads` out to a FASTQ formatted file with the given `filename`.

If this method is used with a paired-end read type, then the FASTQ file will
be interleaved; all R1 reads will be odd records, and all R2 reads will be even
records in the file.

!!! note
    Reads are named according to the sequence in the input genome they came
    from. e.g. `@Reference_1_R1` means the first sequence in the genome, and
    `@Reference_2_R1` means the second sequence in the genome.
"""
function generate(filename::String, reads::Reads)
    open(FASTQ.Writer, filename) do wtr
        generate(wtr, reads)
    end
end

"""
    generate(R1name::String, R2name::String, reads::Reads{<:PairedReads})
    
This method only works for paired reads. Instead of interleaving R1 and R2 reads
in a single FASTQ file, R1 and R2 reads are partitioned into two seperate FASTQ
files.
"""
function generate(R1name::String, R2name::String, reads::Reads{<:PairedReads})
    R1W = open(FASTQ.Writer, R1)
    R2W = open(FASTQ.Writer, R2)
    try
        generate(R1W, R2W, reads)
    finally
        close(R1W)
        close(R2W)
    end
end
generate(wtr::FASTQ.Writer, reads::Reads{<:PairedReads}) = generate(wtr, wtr, reads)


# Read type specializations
# =========================

# Paired end reads
# ----------------

"""
    make_reads(::Type{PairedEnd}, p::MoleculePool, flen::Int, rlen::Int = flen)

Create a set of paired-end reads from a pool of DNA molecules `p`.

`flen` sets the length of forward read, and `rlen` sets the length of the
reverse read. If you only provide `flen`, then the function sets `rlen = flen`.

!!! note
    If a molecule in the pool is not long enough to create a forward and/or
    reverse read, then that molecule will simply be skipped. 
"""
function make_reads(::Type{PairedEnd}, p::MoleculePool, flen::Int, rlen::Int = flen)
    vs = views(p)
    vse = take_paired_ends(vs, flen, rlen)
    return Reads{PairedEnd}(PairedEnd(flen, rlen), genome(p), vse, Dict{Int, Vector{Int}}())
end

summary(x::Reads{PairedEnd}) = string(nreads(x), " paired-end reads:")
function printlen(x::Reads{PairedEnd})
    vs = views(x)
    lf = length(vs[1])
    lr = length(vs[2])
    if lf == lr
        return string(" Read length: ", lf)
    else
        return string(" R1 length: ", lf, ", R2 length: ", lr)
    end
end

function generate(R1W::FASTQ.Writer, R2W::FASTQ.Writer, reads::Reads{PairedEnd})
    vs = views(reads)
    g = genome(reads)
    errs = errors(reads)
    n = 0
    for i in eachindex(vs)
        v = vs[i]
        e = get(errs, i, Int[])
        ref = g[seqid(v)]
        # Extract subsequence from the reference and add errors.
        seq = extract_sequence(ref, v)
        add_errors!(seq, e)
        # Make the quality string for the read, (currently no real meaning!).
        qual = Vector{Int}(undef, length(seq))
        FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, fill(30, length(seq)), qual, 1, length(seq))
        # Create the name for the read...
        fragname = string("Refseq_", seqid(v))
        if isodd(i)
            fragname = string(fragname, "_R1")
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R1W, fqread)
        else
            fragname = string(fragname, "_R2")
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R2W, fqread)
        end
        n += 1
    end
    @info string("- ✔ Wrote ", n, " paired end reads to FASTQ file")
end


# Single end reads
# ----------------

"""
    make_reads(::Type{SingleEnd}, p::MoleculePool, len::Int)

Create a set of single-end reads from a pool of DNA molecules `p`.

`len` sets the length of the reads.

The end (strand) from which the reading begins for each DNA molecule in the
pool is determined at random for each molecule, with 50:50 probability.

If you don't provide a value for `len`, then the function will read each
DNA molecule in it's entirety.

!!! note
    If a molecule in the pool is not long enough to create a forward and/or
    reverse read, then that molecule will simply be skipped. 
"""
function make_reads(::Type{SingleEnd}, p::MoleculePool, len::Int)
    vs = views(p)
    vse = take_single_ends(vs, len)
    return Reads(SingleEnd(len), genome(p), vse, Dict{Int, Vector{Int}}())
end
function make_reads(::Type{SingleEnd}, p::MoleculePool)
    vs = views(p)
    vse = take_single_ends(vs)
    return Reads(SingleEnd(nothing), genome(p), vse, Dict{Int, Vector{Int}}())
end
make_reads(::Type{SingleEnd}, p::MoleculePool, len::Nothing) = make_reads(SingleEnd, p)

summary(x::Reads{<:SingleEnd}) = string(nreads(x), " single-end reads:")
function printlen(x::Reads{SingleEnd{Nothing}})
    mx, av, mn = summarize_lengths(views(x))
    return string(" Maximum read size: $mx\n Average read size: $av\n Minimum read size: $mn")
end
printlen(x::Reads{SingleEnd{UInt}}) = string(" Read length: ", x.tech.len)

total_bp(x::Reads{SingleEnd{UInt}}) = nreads(x) * x.tech.len
total_bp(x::Reads{SingleEnd{Nothing}}) = sum([length(v) for v in views(x)])

function pick_error_positions!(reads::Reads{SingleEnd{UInt}}, nerrs::Int)
    errs = errors(reads)
    empty!(errs)
    len = reads.tech.len
    nr = pick_read(reads)
    for i in 1:nerrs
        readidx = pick_read(reads)
        pos = rand(1:len)
        if haskey(errs, readidx)
            push!(errs[readidx], pos)
        else
            errs[readidx] = [pos]
        end
    end
end

function pick_error_positions!(reads::Reads{SingleEnd{Nothing}}, nerrs::Int)
    errs = errors(reads)
    empty!(errs)
    vs = views(reads)
    for i in 1:nerrs
        readidx = pick_read(reads)
        pos = rand(1:length(vs[readidx]))
        if haskey(errs, readidx)
            push!(errs[readidx], pos)
        else
            errs[readidx] = [pos]
        end
    end
end

function generate(wtr::FASTQ.Writer, reads::Reads{<:SingleEnd})
    vs = views(reads)
    g = genome(reads)
    errs = errors(reads)
    n = 0
    for i in eachindex(vs)
        v = vs[i]
        e = get(errs, i, Int[])
        ref = g[seqid(v)]
        # Extract subsequence from the reference and add errors.
        seq = extract_sequence(ref, v)
        add_errors!(seq, e)
        # Make the quality string for the read, (currently no real meaning!).
        qual = Vector{Int}(undef, length(seq))
        FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, fill(30, length(seq)), qual, 1, length(seq))
        # Create the name for the read...
        fragname = string("Refseq_", seqid(v))
        fqread = FASTQ.Record(fragname, seq, qual)
        write(wtr, fqread)
        n += 1
    end
    @info string("- ✔ Wrote ", n, " single end reads to FASTQ file")
end


# 10x style tagged reads
# ----------------------

"""
    make_reads(::Type{TaggedPairs}, p::MoleculePool, flen::Int, rlen::Int = flen)

Create a set of tagged paired-end reads from a pool of DNA molecules `p`.

`flen` sets the length of forward read, and `rlen` sets the length of the
reverse read. If you only provide `flen`, then the function sets `rlen = flen`.

When a set of `TaggedPairs` is written to file, the tag information is contained
in the R1 read of each read-pair. The first 16bp of each R1 read is a sequence
that is the tag, and a following 7bp are a buffer between the 16bp tag, and the
rest of the read sequence.

!!! note
    If a molecule in the pool is not long enough to create a forward and/or
    reverse read, then that molecule will simply be skipped. 
"""
function make_reads(::Type{TaggedPairs}, p::MoleculePool, flen::Int, rlen::Int = flen)
    tech = TaggedPairs(flen, rlen)
    vs = views(p)
    vse = take_paired_ends(vs, flen - 23, rlen) # 23 is the number of bp of the tag (16), and the buffer (7).
    return Reads{TaggedPairs}(tech, genome(p), vse, Dict{Int, Vector{Int}}())
end

summary(x::Reads{TaggedPairs}) = string(nreads(x), " tagged paired-end reads:")
function printlen(x::Reads{TaggedPairs})
    vs = views(x)
    lf = length(vs[1])
    lr = length(vs[2])
    lf += 23 # 23 is the number of bp of the tag (16), and the buffer (7).
    if lf == lr
        return string(" Read length: ", lf)
    else
        return string(" R1 length: ", lf, ", R2 length: ", lr)
    end
end

function prepare_tags(reads::Reads{TaggedPairs})
    tags = keys(summarize_tags(views(reads)))
    tagrange = 0x0000000000000000:0x00000000ffffffff
    tagseqs = DNAKmer{16}.(sample_values(tagrange, length(tags)))
    tagdict = Dict(zip(tags, tagseqs))
    return tagdict
end

openfastq(f::String) = open(FASTQ.Reader(f))

function generate(R1W::FASTQ.Writer, R2W::FASTQ.Writer, reads::Reads{TaggedPairs})
    
    vs = views(reads)
    g = genome(reads)
    errs = errors(reads)
    
    bufseq = LongSequence{DNAAlphabet{2}}(rand(ACGT, 7))
    tagdict = prepare_tags(reads)
    n = 0
    
    for i in eachindex(vs)
        v = vs[i]
        e = get(errs, i, Int[])
        ref = g[seqid(v)]
        # Extract subsequence from the reference and add errors.
        seq = extract_sequence(ref, v)
        
        if isodd(i) # We are dealing with an R1 read
            tagseq = LongSequence{DNAAlphabet{2}}(tagdict[tag(v)])
            seq = tagseq * bufseq * seq
        end
        
        add_errors!(seq, e)
        # Make the quality string for the read, (currently no real meaning!).
        qual = Vector{Int}(undef, length(seq))
        FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, fill(30, length(seq)), qual, 1, length(seq))
        
        # Create the name for the read...
        fragname = string("Refseq_", seqid(v))
        
        if isodd(i)
            fragname = string(fragname, "_R1")
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R1W, fqread)
        else
            fragname = string(fragname, "_R2")
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R2W, fqread)
        end
        n += 1
    end
    @info string("- ✔ Wrote ", n, " tagged paired end reads to FASTQ file")
end