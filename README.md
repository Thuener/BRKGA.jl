# BRKGA.jl
An implementation in Julia of the Biased Random-Key Genetic Algorithm(http://mauricio.resende.info/doc/srkga.pdf)

To use this package for a specific problem you have to construct a type that is subtype(<:) of BRKGAProblem and implement functions decoder and nalleles. An example of how to create this functions is presented in [TSP.jl](./src/TSP.jl):

```julia
type TSPProblem <: BRKGAProblem
    ...
end

BRKGA.nalleles(tsp::TSPProblem) = npairs(tsp)

function BRKGA.decoder(tsp::TSPProblem, chromosome::Vector{Float64})
    ... 
end
```
