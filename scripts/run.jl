include(joinpath(@__DIR__, "..", "src", "reputation_and_learning.jl"))

# Generate the model instance with the basic parameters
m = SimpleModel()
@info m  # prints the parameters.

# Solve the model
sol = construct_solution(m)

# Print T and m as in the paper
@info "  T : "  sol.T
@info "  m : "  get_m(sol)

# Do Figure 1 in the paper and save it to PDF
p = plots(sol, plot_tmax=80)
savefig(p, joinpath(@__DIR__, "..", "output", "figure.pdf"))

