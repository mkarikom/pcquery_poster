function runLattice(tspan, p, n_vars, dist)
    n_cells = size(dist)[1]
    var_inds = collect(1:n_vars*n_cells)

    u0 = rand(n_cells*n_vars).*10
    params = Dict("p" => p, "dist" => dist, "n_cells" => n_cells, "var_inds" => var_inds, "n_vars" => n_vars)

    function lateralI(du,u,params,t)
        d_ind = reshape(params["var_inds"],2,:)[1,:] # indices of the 1st variable eg delta
        h_ind = reshape(params["var_inds"],2,:)[2,:] # indices of the 2nd variable eg her1
        d_bar = params["dist"]*(u[d_ind]./(ceil.(params["dist"])*ones(params["n_cells"]))) # sum neighbor delta expression
        p = params["p"]
        n_cells = params["n_cells"]
        for ddu = 1:n_cells
            nu = p[(ddu-1)*5+1]
            betaD = p[(ddu-1)*5+2]
            betaR = p[(ddu-1)*5+3]
            h = p[(ddu-1)*5+4]
            m = p[(ddu-1)*5+5]

            du[(ddu-1) * 2 + 1] = nu .* (betaD .* 1 ./ (1 .+ u[h_ind][ddu] .^ h) .- u[d_ind][ddu]) # the value of delta
            du[(ddu-1) * 2 + 2] = (params["dist"]*(betaR .* u[d_ind] .^ m ./ (1 .+ u[d_ind] .^ m)) .- u[h_ind])[ddu] # the value of hes1
        end
    end

    prob = ODEProblem(lateralI,u0,tspan,params)
    #sol = solve(prob,BS3(),tstops = collect(tspan[1]:tspan[2])) # save specific time points, to debug run: "Juno.@enter solve(prob)"
    sol = solve(prob,BS3()) # save specific time points, to debug run: "Juno.@enter solve(prob)"

    plot(sol)

    #=
    Organize the solutions at specified time points, collect the derivatives
    =#

    # reshape the solution output as an array u for each species
    allu = sol.u[1] # |time| x |species| species alternates
    for uu = 2:(size(sol.u)[1])
    	allu = vcat(allu,sol.u[uu])
    end

    u_d = reshape(allu, n_vars, n_cells, size(sol.u)[1])[1,:,:]'
    u_h1 = reshape(allu, n_vars, n_cells, size(sol.u)[1])[2,:,:]'

    #t_ind = findall(x -> x >= 0, sol.t - ceil.(sol.t)) # timepoints where derivatives were automatically saved
    alld = sol(sol.t[1],Val{1})
    for tt = 2:length(sol.t)
    	alld = vcat(alld,sol(sol.t[tt],Val{1})) # get the 1st derivative at time sol.t[tt]
    end
    du_d = reshape(alld, n_vars, n_cells, size(sol.u)[1])[1,:,:]'
    du_h1 = reshape(alld, n_vars, n_cells, size(sol.u)[1])[2,:,:]'

    return u_d, u_h1, du_d, du_h1, sol
end
