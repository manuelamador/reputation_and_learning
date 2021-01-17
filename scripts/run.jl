include(joinpath(@__DIR__, "..", "src", "reputation_and_learning.jl"))

# Generate the model instance with the basic parameters
m = SimpleModel()
@info m  # prints the parameters.

        # Info: SimpleModel{Float64}
        # │   r: Float64 0.15
        # │   i: Float64 0.01
        # │   λ: Float64 0.2
        # │   δ: Float64 0.02
        # │   ϵ: Float64 0.01
        # │   y: Float64 1.0
        # │   bmax: Float64 1.0
        # │   qupperbar: Float64 0.9130434782608696
        # │   qlowerbar: Float64 0.6000000000000001
        # └   cbar: Float64 1.14

# Solve the model
sol = construct_solution(m)

# Print T and m as in the paper
@info "  T : "  sol.T
        # sol.T = 31.17531109613466

@info "  m : "  get_m(sol)
        # get_m(sol) = 197.87246177564631

# Do Figure 1 in the paper and save it to PDF
p = plots(sol, plot_tmax=80)
savefig(p, joinpath(@__DIR__, "..", "output", "figure.pdf"))

