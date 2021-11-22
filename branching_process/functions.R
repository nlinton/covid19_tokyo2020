
# Reporting delay
get_reporting_delay <- \(dat, pars00) {
    rd <- dat$rd
    rd_t <- dat$rd_t
    lik <- \(pars) { return(-sum(log((pweibull(rd+0.5, shape=pars[1], scale=pars[2]) -
                                      pweibull(rd-0.5, shape=pars[1], scale=pars[2])) /
                                      (pweibull(rd_t, shape=pars[1], scale=pars[2]))))) }
    options(warn=-1)
    opt_rd <- optim(pars00, fn=lik, method='BFGS', control = list(maxit=1e8, trace=0))
    options(warn=0)
    opt_rd$par
}

# Backprojection
backproj_once <- \(df, pars1, dists, grp_name, pmf_max=30) {

    T <- nrow(df)

    dist_pmf <- eval(parse(text=paste0(dists, "(1:T, pars1)")))[1:pmf_max]

    # Use surveillance package
    smth_par <- 2
    bp_ctrl <- list(k=smth_par, eps=rep(1e-4,smth_par), iter.max=rep(1000,smth_par), Tmark=nrow(df), B=-1, alpha=0.01, verbose=FALSE, lambda0=NULL, eq3a.method="C")
    sts_out <- surveillance::backprojNP(new("sts", epoch=df$t, observed=df[grp_name]), incu.pmf=dist_pmf, control=bp_ctrl)
    out <- surveillance::upperbound(sts_out); out[out<0.01] <- 0
    df_out <- eval(parse(text=paste0("cbind(df, ",grp_name,"_bp=out)")))

    # Normalize
    df_out |> mutate(!!paste0(grp_name,"_norm") := eval(parse(text=paste0(grp_name,"_bp/sum(",grp_name,"_bp)*sum(",grp_name,")"))))
}

# Create next-generation matrix using a data frame with variables type, geG, G, z0
create_ngm <- \(df, G, D, Rb, theta, grp_names) {
    G <- G # Number of generations
    D <- D # Number of groups
    N_tot <- df %>% group_by(gen) %>% summarize(tot=sum(n)) %>% dplyr::select(tot) %>% pull() # get totals as vector
    N_mat <- matrix(df$n, nrow=G) # N x D matrix of total population by type for each generation
    n_mat <- t(t(N_mat) / N_tot) # divide each cell by corresponding total population N for each generation
    k_by_type <- Rb/(((1-theta) * n_mat) + theta) # solve for k using baseline Rt
    nk_o <- sapply(1:dim(k_by_type)[2], \(i) sapply(1:dim(k_by_type)[1], \(g) k_by_type[g,i]*(1-theta)*n_mat[g,i])) # use k to obtain NGM value for each i in each generation
    ngm_nk_o <- split(t(nk_o), rep(1:ncol(t(nk_o)), each = nrow(t(nk_o)))) # split k for each i into generations
    names(ngm_nk_o) <- paste0("gen", 0:(dim(N_mat)[1]-1)) # name the generations
    ngm_base <- lapply(ngm_nk_o, \(m) matrix(rep(m,D), ncol=D, dimnames=list(grp_names, grp_names))) # change to matrix according to i
    ngm <- lapply(ngm_base, \(g) { diag(g) <- Rb; g } ) # overwrite the diagonal
    return(ngm)
}
