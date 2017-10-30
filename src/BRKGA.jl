# Implementation of the Genetic algorithm BRKGA in Julia
module BRKGA

export Parameters, BRKGAProblem
export start

abstract BRKGAProblem

param(prob::BRKGAProblem) = prob.param
nalleles(prob::BRKGAProblem)  = error("Please implement nalleles function for $(typeof(x)) subtype of BRKGAProblem")
decoder(prob::BRKGAProblem)  = error("Please implement decoder function for $(typeof(x)) subtype of BRKGAProblem")

type Individual
    chromosome::Vector{Float64}
    fitness::Float64
end

function Individual(prob::BRKGAProblem)
    chr = rand(nalleles(prob))
    fit = decoder(prob, chr)
    Individual(chr, fit)
end

function Individual(prob::BRKGAProblem, chr::Vector{Float64})
    fit = decoder(prob, chr)
    Individual(chr, fit)
end

fitness(ind::Individual)          = ind.fitness
chromosome(ind::Individual)       = ind.chromosome
allele(ind::Individual, idx::Int) = ind.chromosome[idx]

type Parameters
    proballeles::Float64 # Crossover probability
    maxpop::Int64        # Max number of the population (has to be bigger than elitepop+ elitepop+ nonelitepop)
    maxgen::Int64        # Max number of generations
    elitepop::Int64      # Number of elite individuals
    randind::Int64       # Number of new random individuals inserted in each generation
    crosstype::Int64     # Crossover type
end

Parameters(proballeles::Float64, maxpop::Int64, maxgen::Int64, felite::Float64, frand::Float64, crosstype::Int64) =
 Parameters(proballeles, maxpop, maxgen, floor(Int64, maxpop*felite), floor(Int64, maxpop*frand), crosstype)

# Util function for parameters
probal(para::Parameters)    = para.proballeles
maxpop(para::Parameters)    = para.maxpop
maxgen(para::Parameters)    = para.maxgen
neli(para::Parameters)      = para.elitepop
crossind(para::Parameters)  = maxpop(para) -(neli(para) + randind(para))
randind(para::Parameters)   = para.randind
crosstype(para::Parameters) = para.crosstype

"""

Make new population crossing-over an elite individual with non-elite to construct a new individual
"""
function crosspop(prob::BRKGAProblem, pop::Vector{Individual})
    para = param(prob)
    newpop = Array(Individual, crossind(para))
    for i = 1:crossind(para)
        # Choose parent
        indpeli = floor(Int64, rand()*neli(para)) +1

        # Choose the other parent
        if crosstype(para) == 1
            # Join elite with non-elite or random
            indpneli = neli(para) + floor(Int64, rand()*(crossind(para)+randind(para))) +1
        elseif crosstype(para) == 2
            # Join elite with elite, non-elite or random
            indpneli = floor(Int64, rand()*maxpop(para)) +1
            # Test if they are the same
            if indpneli == peli
                if indpeli < neli(para)
                    indpeli += 1
                else
                    indpeli -= 1
                end
            end
        end

        peli = pop[indpeli]   # Elite parent
        pneli = pop[indpneli] # Non-elite parent

        # The first parent is always the one with less fitness
        if fitness(peli) > fitness(pneli)
            aux = peli
            peli = pneli
            pneli = aux
        end

        # Joint the parents
        chr = Array(Float64, nalleles(prob))
        for j = 1:nalleles(prob)
            r = rand()

            # Allele form the non-elite parent
            if r > probal(para)
                chr[j] = allele(pneli, j)
            else # Allele form the elite parent
                chr[j] = allele(peli, j)
            end
        end
        newpop[i] = Individual(prob, chr)
    end
    return newpop
end

"""

Create a random population
"""
function randpop(prob::BRKGAProblem, n::Int)
    para = param(prob)
    newpop = Array(Individual, n)
    # Generates a random solution
    for i = 1:n
        newpop[i] = Individual(prob)
    end
    return newpop
end

"""

Start the GA algorithm
"""
function start(prob::BRKGAProblem; targetcost=0, maxstableit=Inf)
    tic()
    startt = time()
    para = param(prob)
    pop = randpop(prob, maxpop(para))

    lastsol = Inf
    g = 0
    stableit = 0 # Iterations without improvement
    cursol = 0
    bestsoltime = 0
    bestsolution = Inf
    for g = 1:maxgen(para)
        sort!(pop, by=fitness)
        cropop = crosspop(prob, pop)
        rndpop = randpop(prob, randind(para))
        pop = vcat(pop[1:neli(para)], cropop, rndpop)
        cursol = fitness(pop[1])

        info("Generation: $g current solution: $(cursol)")

        if cursol < targetcost
            break
        end
        if isapprox(lastsol, cursol)
            info("Iteration without improvement $(stableit)")
            stableit += 1
            if stableit >= maxstableit
                break
            end
        else
            stableit = 0
        end
        if cursol < bestsolution
            bestsoltime = time() -startt
            bestsolution = cursol
        end
        lastsol = cursol
    end
    totaltime = toq()
    info("Final time $(totaltime)")
    return cursol, g, totaltime, bestsoltime
end

end # BRKGA
