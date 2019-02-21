
function sample(vs::Views, n::Int)
    p = resize!(randperm(length(vs)), n)
    return vs[p]
end

needed_sample_size(C::Int, G::Int, L::Int) = div(C * G, L)
expected_coverage(G::Int, L::Int, N::Int) = div(L * N, G)