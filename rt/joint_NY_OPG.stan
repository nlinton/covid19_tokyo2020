functions {
   /*
    * Returns discretized gamma cdf
    *
    * @param alpha  Alpha parameter of gamma distribution
    * @param beta   Beta parameter of gamma distribution
    * @param K      Length of input vector
    *
    * return cdf
    */
    vector pgamma(real alpha, real beta, int K) {
      vector[K] res;
      for (k in 1:K)
        res[k] = gamma_cdf(k, alpha, beta);
        return append_row(res[1], res[2:K]-res[1:(K-1)]);
    }

   /*
    * Returns convolution of two vectors
    *
    * @param X      Vector X
    * @param Yrev   Reverse of vector Y
    * @param K      Length of vectors X and Yrev
    *
    * return cdf
    */
    vector convolution(vector X, vector Yrev, int K) {
        vector[K-1] res;
        res[1] = X[1]*Yrev[K];
        for (k in 2:K-1)
            res[k] = dot_product(head(X, k), tail(Yrev, k));
        return res;
    }

   /*
    * Returns assortativity-based next-generation matrix
    *
    *
    * return cdf
    */
    matrix get_ngm(real theta, vector n, vector Rt, int G, int D) {
        matrix[G, G] res;
        res = [[ Rt[1], ( (1.0 - theta) * n[2] ) * ( Rt[2] / ( theta + (1.0 - theta) * n[2] ) ) ], // first row of ngm
                [ ( (1.0 - theta) * n[1] ) * ( Rt[1] / ( theta + (1.0 - theta) * n[1] ) ),  Rt[2] ]]; // second row of ngm
        return res;
    }
}

data {
   int<lower=1> D_ny; // number of days in NY dataset
   int<lower=1> D_opg; // number of days in Olympic and Paralympic Games (OPG) dataset
   int<lower=1> T_ny; // number of Rt steps for NY time period
   int<lower=1> T_opg; // number of Rt steps for OPG time period
   int<lower=1> G; // number of groups during OPG period (OPG, Tokyo)

   // for assortativity
   real<lower=0, upper=1> theta; // assortativity measure
   vector<lower=0>[G] n; // population size of each group

   // externals of renewal equation
   int<lower=1> vec_ny[T_ny];
   int<lower=1> vec_opg[T_opg];
   int<lower=0> t_ny1[vec_ny[1]]; 
   int<lower=0> t_ny2[vec_ny[2]]; 
   int<lower=0> t_opg1[vec_opg[1]]; 
   int<lower=0> t_opg2[vec_opg[2]]; 
   int<lower=0> t_opg3[vec_opg[3]]; 
   int<lower=0> t_opg4[vec_opg[4]]; 
   int<lower=0> t_opg5[vec_opg[5]]; 

   // Considering spectators
   int<lower=0> t_ny1c[vec_ny[1]]; 
   int<lower=0> t_ny2c[vec_ny[2]]; 
   int<lower=0> t_opg1c[vec_opg[1]]; 
   int<lower=0> t_opg2c[vec_opg[2]];
   int<lower=0> t_opg3c[vec_opg[3]]; 
   int<lower=0> t_opg4c[vec_opg[4]];
   int<lower=0> t_opg5c[vec_opg[5]];
   real<lower=0, upper=1> k_frac_g;
   real<lower=0, upper=1> k_frac_t;

   // internals of renewal equation
   vector<lower=0>[D_ny] cases_ny;
   matrix<lower=0>[D_opg, G] cases_opg;
   real<lower=0> gi_par1;
   real<lower=0> gi_par2;
}

transformed data {
   vector[D_ny-1] conv_ny; 
   array[G] vector[D_opg-1] conv_opg; 
   {
   // Convolution NY
   vector[D_ny] gi_ny = pgamma(gi_par1, gi_par2, D_ny);
   vector[D_ny] gi_rev_ny;
   for (d in 1:D_ny) 
       gi_rev_ny[d] = gi_ny[D_ny+1-d];
   conv_ny = convolution(cases_ny, gi_rev_ny, D_ny); 

   // Convolution OPG
   vector[D_opg] gi_opg = pgamma(gi_par1, gi_par2, D_opg); 
   vector[D_opg] gi_rev_opg;
   for (d in 1:D_opg) 
       gi_rev_opg[d] = gi_opg[D_opg+1-d];
   for (g in 1:G) 
       conv_opg[g] = convolution(cases_opg[,g], gi_rev_opg, D_opg); 
   }
}

parameters {
   real<lower=0> rel_k;
   real<lower=0> Rt_ny;
   array[T_opg] vector<lower=0>[G] Rt_opg;
}

transformed parameters {
    // Original fit
    matrix<lower=0>[G, G] ngm1 = get_ngm(theta, n, Rt_opg[1], G, D_opg);
    matrix<lower=0>[G, G] ngm2 = get_ngm(theta, n, Rt_opg[2], G, D_opg);
    matrix<lower=0>[G, G] ngm3 = get_ngm(theta, n, Rt_opg[3], G, D_opg);
    matrix<lower=0>[G, G] ngm4 = get_ngm(theta, n, Rt_opg[4], G, D_opg);
    matrix<lower=0>[G, G] ngm5 = get_ngm(theta, n, Rt_opg[5], G, D_opg);

    // Spectators counterfactual
    real<lower=0> rel_k_g = ((rel_k - 1) * k_frac_g) + 1; // Games spectators counterfactual rel_k
    real<lower=0> rel_k_t = ((rel_k - 1) * k_frac_t) + 1; // Tokyo spectators counterfactual rel_k
    matrix<lower=0>[G, G] ngm2sg = get_ngm(theta, n, Rt_opg[2] * rel_k_g, G, D_opg); // Spectator counterfactual
    matrix<lower=0>[G, G] ngm3sg = get_ngm(theta, n, Rt_opg[3] * rel_k_g, G, D_opg); // Spectator counterfactual
    matrix<lower=0>[G, G] ngm5sg = get_ngm(theta, n, Rt_opg[5] * rel_k_g, G, D_opg); // Spectator counterfactual
    matrix<lower=0>[G, G] ngm2st = get_ngm(theta, n, Rt_opg[2] * rel_k_t, G, D_opg); // Spectator counterfactual
    matrix<lower=0>[G, G] ngm3st = get_ngm(theta, n, Rt_opg[3] * rel_k_t, G, D_opg); // Spectator counterfactual
    matrix<lower=0>[G, G] ngm5st = get_ngm(theta, n, Rt_opg[5] * rel_k_t, G, D_opg); // Spectator counterfactual
}

model {
    // Fit to NY data
    target += gamma_lpdf(cases_ny[t_ny1] | Rt_ny * conv_ny[t_ny1c], 1.0) + log(1e-8)
            + gamma_lpdf(cases_ny[t_ny2] | Rt_ny * rel_k * conv_ny[t_ny2c], 1.0) + log(1e-8);

    // Fit to OPG data
    for (g in 1:G) {
        target += gamma_lpdf(cases_opg[t_opg1,g] | ngm1[g, 1] * conv_opg[g][t_opg1c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg1,g] | ngm1[g, 2] * conv_opg[g][t_opg1c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg2,g] | ngm2[g, 1] * conv_opg[g][t_opg2c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg2,g] | ngm2[g, 2] * conv_opg[g][t_opg2c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg3,g] | ngm3[g, 1] * conv_opg[g][t_opg3c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg3,g] | ngm3[g, 2] * conv_opg[g][t_opg3c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg4,g] | ngm4[g, 1] * conv_opg[g][t_opg4c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg4,g] | ngm4[g, 2] * conv_opg[g][t_opg4c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg5,g] | ngm5[g, 1] * conv_opg[g][t_opg5c], 1.0) + log(1e-8)
                + gamma_lpdf(cases_opg[t_opg5,g] | ngm5[g, 2] * conv_opg[g][t_opg5c], 1.0) + log(1e-8);
    }
}

generated quantities {

    vector<lower=0>[D_ny-14] Ect_ny = append_row(Rt_ny * conv_ny[t_ny1c], Rt_ny * rel_k * conv_ny[t_ny2c]);

    real<lower=0> Rt_ny_elev = Rt_ny * rel_k;
    array[T_opg] vector<lower=0>[G] Rt_spec = { [ Rt_opg[1,1],           Rt_opg[1,2] ]',
                                                [ Rt_opg[2,1] * rel_k_g, Rt_opg[2,2] * rel_k_t ]',
                                                [ Rt_opg[3,1] * rel_k_g, Rt_opg[3,2] * rel_k_t ]',
                                                [ Rt_opg[4,1],           Rt_opg[4,2] ]',
                                                [ Rt_opg[5,1] * rel_k_g, Rt_opg[5,2] * rel_k_t ]' };

    array[G] vector<lower=0>[D_opg-14] Ect_opg;
    for (g in 1:G) {
        Ect_opg[g] = append_row((ngm1[g, 1] * conv_opg[g][t_opg1c]) + (ngm1[g, 2] * conv_opg[g][t_opg1c]),
                     append_row((ngm2[g, 1] * conv_opg[g][t_opg2c]) + (ngm2[g, 2] * conv_opg[g][t_opg2c]),
                     append_row((ngm3[g, 1] * conv_opg[g][t_opg3c]) + (ngm3[g, 2] * conv_opg[g][t_opg3c]),
                     append_row((ngm4[g, 1] * conv_opg[g][t_opg4c]) + (ngm4[g, 2] * conv_opg[g][t_opg4c]),
                                (ngm5[g, 1] * conv_opg[g][t_opg5c]) + (ngm5[g, 2] * conv_opg[g][t_opg5c]))))); 
    }

    // Use transformed spectator NGM during Games, use fitted NGM during other times
    array[G] vector<lower=0>[D_opg-14] Ect_spec;
        Ect_spec[1] = append_row((ngm1[1, 1] * conv_opg[1][t_opg1c])   + (ngm1[1, 2]   * conv_opg[1][t_opg1c]),
                      append_row((ngm2sg[1, 1] * conv_opg[1][t_opg2c]) + (ngm2sg[1, 2] * conv_opg[1][t_opg2c]),
                      append_row((ngm3sg[1, 1] * conv_opg[1][t_opg3c]) + (ngm3sg[1, 2] * conv_opg[1][t_opg3c]),
                      append_row((ngm4[1, 1] * conv_opg[1][t_opg4c])   + (ngm4[1, 2]   * conv_opg[1][t_opg4c]),
                                 (ngm5sg[1, 1] * conv_opg[1][t_opg5c]) + (ngm5sg[1, 2] * conv_opg[1][t_opg5c])))));
        Ect_spec[2] = append_row((ngm1[2, 1] * conv_opg[2][t_opg1c])   + (ngm1[2, 2]   * conv_opg[2][t_opg1c]),
                      append_row((ngm2st[2, 1] * conv_opg[2][t_opg2c]) + (ngm2st[2, 2] * conv_opg[2][t_opg2c]),
                      append_row((ngm3st[2, 1] * conv_opg[2][t_opg3c]) + (ngm3st[2, 2] * conv_opg[2][t_opg3c]),
                      append_row((ngm4[2, 1] * conv_opg[2][t_opg4c])   + (ngm4[2, 2]   * conv_opg[2][t_opg4c]),
                                 (ngm5st[2, 1] * conv_opg[2][t_opg5c]) + (ngm5st[2, 2] * conv_opg[2][t_opg5c]))))); 
}
