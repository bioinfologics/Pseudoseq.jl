"Create an empty blueprint for `n` copies of a chromosome of `len` base pairs."
plan_chrom(len::Int, n::Int) = ChromosomeBlueprint(len, n)

function plan_chrom(; len::Int = 500, ploidy::Int = 1, ngenome::Int = 1, het = 0.0, div::Float64 = 0.0)
    bp = ChromosomeBlueprint(len, ploidy * ngenome)
    @info "- ✔ Created blank blueprint"
    if div > 0.0
        if ngenome > 1
            nsites = Int(ceil(len * div))
            bp = plan_het(bp, nsites, repeat(1:ngenome, inner = ploidy))
            @info "- ✔ Planned genome diverging polymorphisms"
        else
            @error "- ✘ Could not plan genome diverging polymorphisms: only one genome is planned"
        end
    end
    
    bp = _try_plan_het(bp, ploidy, ngenome, het)
    
    return bp
end

function _try_plan_het(bp::ChromosomeBlueprint, ploidy::Int, ngenome::Int, het::Float64)
    @info "A single value $het has been provided to plan heterozygosity"
    if het > 0.0 && ngenome === 1 && ploidy === 2
        nsites = Int(floor(length(bp) * het))
        bp = plan_het(bp, nsites, 2)
        @info "- ✔ Planned heterozygosity"
    else
        @error """- ✘ Could not plan heterozygosity:
                      Don't know how to plan heterozygosity when:
                      Ploidy is $ploidy
                      The number of genomes is $ngenome
                      het has been set as $het"""
    end
    return bp
end