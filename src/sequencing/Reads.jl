# Paired end reads
# ================

abstract type Pairing end
struct Paired <: Pairing end
struct UnPaired <: Pairing end

struct Substitution
    pos::UInt64
    base::DNA
end

struct Reads{P<:Pairing,V<:AbstractSequencingView}
    genome::Vector{LongSequence{DNAAlphabet{2}}}
    views::Vector{V}
    #errors::Dict{Int,Vector{Int}}
    substitutions::Dict{Int,Vector{Substitution}}
end

@inline views(x::Reads) = x.views
@inline genome(x::Reads) = x.genome
#@inline errors(x::Reads) = x.errors
@inline substitutions(x::Reads) = x.substitutions
@inline nreads(x::Reads) = length(x.views)
@inline pick_read(reads::Reads) = rand(1:nreads(reads))
#@inline nerrors(x::Reads) = length(x.errors) > 0 ? sum(length(x) for x in values(x.errors)) : 0
@inline function nsubstitutions(x::Reads)
    subs = substitutions(x)
    return isempty(subs) ? 0 : sum(length(x) for x in values(subs))
end

# I/O and printing functions
# --------------------------

summary(x::Reads{Paired,T}) where {T} = string(nreads(x), " paired reads:")
summary(x::Reads{UnPaired,T}) where {T} = string(nreads(x), " unpaired reads:")

function printlen(x::Reads)
    mx, av, mn = summarize_lengths(views(x))
    return string(" Maximum read size: $mx\n Average read size: $av\n Minimum read size: $mn")
end

function Base.show(io::IO, reads::Reads)
    println(summary(reads))
    println(printlen(reads))
    #println(string(" Number of errors: ", nerrors(reads)))
    println(string(" Number of errors: ", nsubstitutions(reads)))
end

# Constructors
# ------------

function Reads{Paired,V}(p::Molecules{V}, flen::Int, rlen::Int = flen) where {V<:AbstractSequencingView}
    vs = views(p)
    vse = take_paired_ends(vs, flen, rlen)
    rds = Reads{Paired,V}(genome(p), vse, Dict{Int,Vector{Substitution}}()) 
    return rds
end

function Reads{UnPaired,V}(p::Molecules{V}, len::Union{Int,Nothing}) where {V<:AbstractSequencingView}
    vs = views(p)
    vse = take_single_ends(vs, len)
    rds = Reads{UnPaired,V}(genome(p), vse, Dict{Int, Vector{Substitution}}())
    return rds
end

"""
    paired_reads(p::MoleculePool, flen::Int, rlen::Int = flen)

Create a set of paired-end reads from a pool of DNA molecules `p`.

`flen` sets the length of forward read, and `rlen` sets the length of the
reverse read. If you only provide `flen`, then the function sets `rlen = flen`.

!!! note
    If a molecule in the pool is not long enough to create a forward and/or
    reverse read, then that molecule will simply be skipped. 
"""
function paired_reads(p::Molecules{T}, flen::Int, rlen::Int = flen) where {T<:AbstractSequencingView}
    return Reads{Paired,T}(p, flen, rlen)
end

"""
    unpaired_reads(p::Molecules{T}, len::Int) where {T<:AbstractSequencingView}

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
function unpaired_reads(p::Molecules{T}, len::Union{Int,Nothing}) where {T<:AbstractSequencingView}
    return Reads{UnPaired,T}(p, len)
end

#=
struct LinearDecayErrorRate
    start::Float64
    decrement::Float64
end

function (x::LinearDecayErrorRate)(rates::Vector{Float64}, sequence::LongSequence{DNAAlphabet{2}})
    for i in eachindex(rates)
        rates[i] = x.start - (i - 1) * x.decrement
    end
end
=#

const posbases = Dict{DNA, Tuple{DNA,DNA,DNA}}(
    DNA_A => (DNA_C, DNA_G, DNA_T),
    DNA_C => (DNA_A, DNA_G, DNA_T),
    DNA_G => (DNA_A, DNA_C, DNA_T),
    DNA_T => (DNA_A, DNA_C, DNA_G)
)

const posbases2 = (
    (DNA_C, DNA_G, DNA_T),
    (DNA_A, DNA_G, DNA_T),
    (DNA_N, DNA_N, DNA_N),
    (DNA_A, DNA_C, DNA_T),
    (DNA_N, DNA_N, DNA_N),
    (DNA_N, DNA_N, DNA_N),
    (DNA_N, DNA_N, DNA_N),
    (DNA_A, DNA_C, DNA_G)
)

mutable struct RandomSubstitutionScatter <: Function
    errs_per_read::Dict{Int, Int}
    readidx::Int
end

function RandomSubstitutionScatter(nsubs::Int, nreads::Int)
    read_idxs = rand(Base.OneTo(nreads), nsubs)
    errs_per_read = Dict{Int, Int}()
    for idx in read_idxs
        errs_per_read[idx] = get(errs_per_read, idx, 0) + 1
    end
    return RandomSubstitutionScatter(errs_per_read, 0)
end

function (f::RandomSubstitutionScatter)(output::Vector{Substitution}, readseq::LongSequence{DNAAlphabet{2}})
    current_read = f.readidx + 1
    f.readidx = current_read
    nerrors = get(f.errs_per_read, current_read, 0)
    posqueue = Random.shuffle(Base.OneTo(length(readseq)))
    resize!(output, nerrors)
    for i in 1:nerrors
        basepos = posqueue[i]
        basenuc = readseq[basepos]
        output[i] = Substitution(basepos, rand(posbases[basenuc])) 
    end
end

struct ClearSubstitutions <: Function end

function (f::ClearSubstitutions)(output::Vector{Substitution}, readseq::LongSequence{DNAAlphabet{2}})
    return empty!(output)
end

struct FixedProbSubstitutions <: Function
    prob::Float64
end

function (f::FixedProbSubstitutions)(output::Vector{Substitution}, readseq::LongSequence{DNAAlphabet{2}})
    p = f.prob
    dice = rand(length(readseq))
    dosub = dice .< p
    @inbounds for pos in 1:length(readseq)
        if dosub[pos]
            basenuc = reinterpret(UInt8, readseq[pos])
            push!(output, Substitution(pos, rand(posbases2[basenuc])))
        end
    end
end

function edit_substitutions!(f::Function, reads::Reads)
    vs = views(reads)
    subs = substitutions(reads)
    subsbuf = Vector{Substitution}()
    for (i, v) in enumerate(vs)
        readseq = extract_sequence(genome(reads), v)
        readsubs = get(subs, i, subsbuf)
        f(readsubs, readseq)
        if isempty(readsubs)           # Read no longer has any subs.
            delete!(subs, i)
        elseif readsubs === subsbuf    # Read did not have subs but now does.
            subs[i] = copy(subsbuf)
            empty!(subsbuf)
        else                           # Read did have subs, and still does.
            subs[i] = readsubs
        end
    end
end

#=
function pick_error_positions!(reads::Reads, nerrs::Int)
    errs = errors(reads)
    empty!(errs)
    vs = views(reads)
    for i in 1:nerrs
        readidx = rand(1:nreads(reads))
        read = vs[readidx]
        pos = rand(1:length(read))
        if haskey(errs, readidx)
            push!(errs[readidx], pos)
        else
            errs[readidx] = [pos]
        end
    end
end

function mark_errors!(reads::Reads, rate::Float64)
    nbases = sum(length.(views(reads)))
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
    newreads = typeof(reads)(genome(reads), views(reads), Dict{Int, Vector{Int}}())
    mark_errors!(newreads, rate)
    return newreads
end
=#

# Read file generation
# --------------------

function add_errors!(seq::BioSequence, errs::Vector{Int})
    for err in errs
        truebase = seq[err]
        seq[err] = rand(posbases[truebase])
    end
end

function generate(wtr::FASTQ.Writer, reads::Reads{UnPaired,BasicSequencingView})
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

function generate(R1W::FASTQ.Writer, R2W::FASTQ.Writer, reads::Reads{Paired,BasicSequencingView})
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
        if isodd(i)
            fragname = string("readpair_", i)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R1W, fqread)
        else
            fragname = string("readpair_", i - 1)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R2W, fqread)
        end
        n += 1
    end
    @info string("- ✔ Wrote ", n, " paired end reads to FASTQ file")
end

function prepare_tags(reads::Reads{Paired,TaggedSequencingView})
    tags = keys(summarize_tags(views(reads)))
    tagrange = 0x0000000000000000:0x00000000ffffffff
    tagseqs = DNAMer{16}.(sample_values(tagrange, length(tags)))
    tagdict = Dict(zip(tags, tagseqs))
    return tagdict
end

function generate(R1W::FASTQ.Writer, R2W::FASTQ.Writer, reads::Reads{Paired,TaggedSequencingView})
    
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
        #fragname = string("Refseq_", seqid(v))
        
        if isodd(i)
            fragname = string("readpair_", i)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R1W, fqread)
        else
            fragname = string("readpair_", i - 1)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R2W, fqread)
        end
        n += 1
    end
    @info string("- ✔ Wrote ", n, " tagged paired end reads to FASTQ file")
end

"""
    generate(R1name::String, R2name::String, reads::Reads{Paired,<:AbstractSequencingView})
    
This method only works for paired reads. Instead of interleaving R1 and R2 reads
in a single FASTQ file, R1 and R2 reads are partitioned into two seperate FASTQ
files.
"""
function generate(R1name::String, R2name::String, reads::Reads{Paired,<:AbstractSequencingView})
    R1W = open(FASTQ.Writer, R1name)
    R2W = open(FASTQ.Writer, R2name)
    try
        generate(R1W, R2W, reads)
    finally
        close(R1W)
        close(R2W)
    end
end

generate(wtr::FASTQ.Writer, reads::Reads{Paired,<:AbstractSequencingView}) = generate(wtr, wtr, reads)

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






#=

# Read type specializations
# =========================

# Paired end reads
# ----------------



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
        if isodd(i)
            fragname = string("readpair_", i)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R1W, fqread)
        else
            fragname = string("readpair_", i - 1)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R2W, fqread)
        end
        n += 1
    end
    @info string("- ✔ Wrote ", n, " paired end reads to FASTQ file")
end


# Single end reads
# ----------------
#=
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
=#
function make_reads(::Type{SingleEnd}, p::MoleculePool)
    vs = views(p)
    vse = take_single_ends(vs)
    return Reads(SingleEnd(nothing), genome(p), vse, Dict{Int, Vector{Int}}())
end
make_reads(::Type{SingleEnd}, p::MoleculePool, len::Nothing) = make_reads(SingleEnd, p)



total_bp(x::Reads{SingleEnd{UInt}}) = nreads(x) * x.tech.len
total_bp(x::Reads{SingleEnd{Nothing}}) = sum([length(v) for v in views(x)])



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
    tagseqs = DNAMer{16}.(sample_values(tagrange, length(tags)))
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
        #fragname = string("Refseq_", seqid(v))
        
        if isodd(i)
            fragname = string("readpair_", i)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R1W, fqread)
        else
            fragname = string("readpair_", i - 1)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R2W, fqread)
        end
        n += 1
    end
    @info string("- ✔ Wrote ", n, " tagged paired end reads to FASTQ file")
end
=#