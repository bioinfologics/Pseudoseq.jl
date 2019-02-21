
function multiply_views(views::Views, n::Integer)
    newviews = Views(length(views) * n)
    i = 0
    for view in views
        @inbounds for _ in 1:n
            i += 1
            newviews[i] = view
        end
    end
    return newviews
end
