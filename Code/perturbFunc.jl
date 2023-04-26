function runPerturb(p, n_cells, trim_ends, mult, myalpha, tspan, n_runs, n_vars, adj, dist)
    U_d, U_h1, dU_d, dU_h1, sol = runLattice(tspan, p, n_vars, dist)
    t = sol.t

    for run = 2:n_runs
        u_d, u_h1, du_d, du_h1, sol = runLattice(tspan, p, n_vars, dist)
        U_d = vcat(U_d, u_d)
        U_h1 = vcat(U_h1, u_h1)
        dU_d = vcat(dU_d, du_d)
        dU_h1 = vcat(dU_h1, du_h1)
        t = vcat(t, sol.t)
    end

    #perturb the data with some gaussian noise
    U_d_orig = U_d
    U_h1_orig = U_h1
    dU_d_orig = dU_d
    dU_h1_orig = dU_h1

    U_d = U_d + U_d .* mult.*randn(size(U_d,1),size(U_d,2))
    U_h1 = U_h1 + U_h1 .* mult.*randn(size(U_h1,1),size(U_h1,2))
    dU_d = dU_d + dU_d .* mult.*randn(size(dU_d,1),size(dU_d,2))
    dU_h1 = dU_h1 + dU_h1 .* mult.*randn(size(dU_h1,1),size(dU_h1,2))


    results = Dict("u" => Dict("Delta" => U_d', "Hes1" => U_h1'), "du" => Dict("Delta" => dU_d', "Hes1" => dU_h1'), "timepoints" => t, "adj" => adj) # these matrices get transposed when loaded in matlab for some reason

    # loop over all the variables and run dictionary elastic net

    CV =  Array{Float64}(undef, 2 * n_cells + 1, n_cells)
    for nn = 1:n_cells
        h3_dict = U_d[1:end-trim_ends,:] .^ 3 ./ (1 .+ U_d[1:end-trim_ends,:] .^ 3)
        h4_dict = U_d[1:end-trim_ends,:] .^ 4 ./ (1 .+ U_d[1:end-trim_ends,:] .^ 4)
        X = dict = hcat(U_h1[1:end-trim_ends,nn],h3_dict,h4_dict)
        y = dU_h1[1:end-trim_ends,nn]
        cv = glmnetcv(X, y, alpha = myalpha)
        lam = argmin(cv.meanloss)
        CV[:,nn] = cv.path.betas[:, lam]
    end

    CV_orig =  Array{Float64}(undef, n_cells + 1, 4*n_cells)
    for nn = 1:n_cells
        # h2_dict = U_d .^ 2 ./ (1 .+ U_d .^ 2)
        # h3_dict = U_d .^ 3 ./ (1 .+ U_d .^ 3)
        # h4_dict = U_d .^ 4 ./ (1 .+ U_d .^ 4)
        #
        # X = dict = hcat(U_h1,h2_dict,h3_dict,h4_dict)
        # y = dU_h1[:,nn]
        # cv = glmnetcv(X, y, alpha = 0.5)
        # global CV_orig[:,nn] = cv.path.betas[:,end]
    end

    ADJ = ceil.(abs.(CV)) #this is the estimated adjacency
    errors = vcat(2*adj,2*adj) .- ADJ[2:end,:]

    fn = length(findall(x->x == 2, errors))
    fp = length(findall(x->x == -1, errors))
    tp = length(findall(x->x == 1, errors))
    tn = length(findall(x->x == 0, errors))

    tpr = tp/(tp+fn)
    tnr = tn/(tn+fp)
    fnr = fn/(fn+tp)
    fpr = fp/(fp+tn)

    acc = (tp + tn) / (tp + tn + fp + fn)
    print([acc;tpr;tnr;fnr;fpr]')

    nt = length(sol.t)
    samp_U_d = U_d[end-nt:end,:]
    samp_U_h1 = U_h1[end-nt:end,:]
    samp_dU_d = dU_d[end-nt:end,:]
    samp_dU_h1 = dU_h1[end-nt:end,:]

    scatter(sol.t, samp_U_d[1:end-1,:], legend = false); xlabel!("time"); ylabel!("Dll1 mRNA")
    savefig(string("/home/au/mkarikom/Lab/My Talks/NIPS/Poster/Image/dll_17cell_alpha-", myalpha, "_gamma-", mult, ".pdf"))
    scatter(sol.t, samp_U_h1[1:end-1,:], legend = false); xlabel!("time"); ylabel!("Hes1 mRNA")
    savefig(string("/home/au/mkarikom/Lab/My Talks/NIPS/Poster/Image/hes1_17cell_alpha-", myalpha, "_gamma-", mult, ".pdf"))
    plot(sol.t, U_h1_orig[end-nt:end-1,:], legend = false); xlabel!("time"); ylabel!("Hes1 mRNA")
    plot(sol.t, U_d_orig[end-nt:end-1,:], legend = false); xlabel!("time"); ylabel!("Dll1 mRNA")
    return([acc, tpr, tnr, fnr, fpr])
end
