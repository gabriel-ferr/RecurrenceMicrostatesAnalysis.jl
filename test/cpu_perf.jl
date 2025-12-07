using RecurrenceMicrostatesAnalysis

function main()
    data = rand(10000, 10000)
    probs = []

    for i in 1:10000
        push!(probs, distribution(StateSpaceSet(data[i, :]), 0.27, 4))
    end

    return probs
end

main()