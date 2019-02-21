function summarize_lengths(vs::Views, nbins::Int = 50)
    lens = [length(v) for v in vs]
    show(histogram(lens, nbins = nbins, closed = :left))
    println(string("\nMaximum fragment length: ", maximum(lens)))
    println(string("Average fragment length: ", div(sum(lens), length(lens))))
    println(string("Minimum fragment length: ", minimum(lens)))
end