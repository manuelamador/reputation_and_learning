# Simulation for Reputation and Sovereign Default 

This file contains the code that generates Figure 1 in the paper "Reputation and 
Sovereign Default". It also prints the value of T and the value of m, as defined 
in the paper. 

## Running the code 

The code is in Julia. From the root of this folder run julia and: 

    using Pkg
    Pkg.activate(".")
    Pkg.instantiate(".")

    include(joinpath("scripts", "run.jl"))

The PDF figure is saved in the `output` folder. 

