
function tag_views(vs::Views, npools::Int)
    newviews = Views(length(vs))
    possibletags = rand(UInt64, npools)
    tags = rand(possibletags, length(vs))
    for i in eachindex(vs)
        view = vs[i]
        newviews[i] = SequencingView(view.id, first(view), last(view), tags[i])
    end
    return newviews
end