function summarize_lengths(vs::Views)::Tuple{Int, Int, Int}
    lens = [length(v) for v in vs]
    return maximum(lens), div(sum(lens), length(vs)), minimum(lens)
end

function summarize_tags(vs::Views)
    tagdict = Dict{TAG_TYPE, Int}()
    for v in vs
        t = tag(v)
        if t > 0
            tagdict[t] = get(tagdict, t, 0) + 1
        end
    end
    return tagdict
end