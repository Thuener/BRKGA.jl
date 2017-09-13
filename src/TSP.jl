# TSP instance for using BRKGA using tsplib
module TSP

using Distances, BRKGA, Logging

export TSPProblem
export readTSP, decoder
export npairs, pair, param, nalleles

type TSPProblem <: BRKGAProblem
    pairs::Vector{Vector{Float64}}
    param::Parameters
end

# Funtions for TSPProblem
npairs(tsp::TSPProblem)         = length(tsp.pairs)
pair(tsp::TSPProblem, idx::Int) = tsp.pairs[idx]
BRKGA.nalleles(tsp::TSPProblem) = npairs(tsp)

"""

Read the TSPlib file
"""
function readTSP(file::String)
    pairs = Vector{Vector{Float64}}
    open(file,"r") do f
        readline(f) # NAME:
        readline(f) # COMMENT:
        readline(f) # TYPE:

        line = readline(f) # DIMENSION:
        items = split(line, " ")
        dim = parse(Int, items[2])
        readline(f) # EDGE_WEIGHT_TYPE
        readline(f) #NODE_COORD_SECTION

        pairs = Array(Vector{Float64}, dim)
        for i = 1:dim
            line = readline(f)
            items = matchall(r"[0-9]+", line)
            x = parse(Int64, items[end-1])
            y = parse(Int64, items[end])
            pairs[i] = [x, y]
        end
    end
    return pairs
end

"""

Evaluate distance between two nodes
"""
function distance(p1::Vector{Float64}, p2::Vector{Float64})
    return floor(Int64, euclidean(p1, p2)+0.5)
end

"""

Decoder for TSP
"""
function BRKGA.decoder(tsp::TSPProblem, chromosome::Vector{Float64})

    idx = sortperm(chromosome)

    # Evaluate the tour distance
    dist = 0
    for i = 2:npairs(tsp)
        p1 = pair(tsp, idx[i-1])
        p2 = pair(tsp, idx[i])
        dist += distance(p1, p2)
    end
    firt = pair(tsp, idx[1])
    last = pair(tsp, idx[end])
    dist += distance(firt, last)
    return dist
end

end #TSP
