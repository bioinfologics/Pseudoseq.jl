struct Universe
    seqviews::Views
end

function makeuniverse(rdr::FASTA.Reader, ng::Int, iscircular::Bool = false)
    @info "Creating universe..."
    vs = Views(rdr, iscircular)
    @info "Populating universe with $ng copies of the genome..."
    mvs = multiply_views(vs, ng)
    @info "Done"
    return Universe(mvs)
end

"""
    makeuniverse(file::String, n::Int, iscircular::Bool = false)

Create a universe of `ng` copies of the genome in the fasta formatted `file`.
"""
function makeuniverse(file::String, ng::Int, iscircular::Bool = false)
    open(FASTA.Reader, file) do rdr
        makeuniverse(rdr, n, iscircular)
    end
end

views(u::Universe) u.seqviews