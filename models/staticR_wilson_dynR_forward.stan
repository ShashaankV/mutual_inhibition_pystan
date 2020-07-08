// fits static rivalry and intermittent (dynamic S)
// assume shared parameters except gamma
// make sure tau are consistent, tau_a in Wilson for static is unit magnitude
// dynamic S has tau_u and tau_a, so tau_u is set to tau_a/1000
// dt is tau_u/10

functions
{
#include ../models/functions/wilson_fun.stan
#include ../models/functions/dynamicS_fun.stan
}

data 
    {
        real J; 
        real beta;
    
        int<lower = 0> Ns; // number of static rivalry measurements
        real c[Ns,2]; //contrast (stimulus) condition
        real x_r[1]; // junk for solver
        int x_i[1]; // junk
        int max_steps;
        real rel_tol;
        real function_tol;

        int<lower = 0> Nd; // number of intermittent rivalry measurements
                             // observation, report noise;
}

transformed data
    {
        real cvd;
        cvd = 0.6;
    }

generated quantities //
    {   
        real Ts_obs[Ns, 2]; // synthetic T static 
        real uinf[Ns,2];
        real t12[1,2];

        real lambda;
        real Tdyn_obs[Nd]; //synthetic noisy report during intermittent rivalry 

        for(n in 1:Ns)
        {
            vector[2] uinf_input_1 = [c[n,1],J]';
            vector[2] uinf_input_2 = [c[n,2],J]';
            vector[1] guess = [0.01]';
            vector[1] sol_1 = algebra_solver(f, guess, uinf_input_1, x_r, x_i,rel_tol,function_tol,max_steps);
            vector[1] sol_2 = algebra_solver(f, guess, uinf_input_2, x_r, x_i,rel_tol,function_tol,max_steps);
            uinf[n,1] = sol_1[1];
            uinf[n,2] = sol_2[1];

            // u1logarg[n] = (J*c[n,1]^0.1 - beta*uinf[n,2])/(uinf[n,2]);
            // u2logarg[n] = (J*c[n,2]^0.1 - beta*uinf[n,1])/(uinf[n,1]);

            t12 = T12fun(c[n,1],c[n,2],uinf[n,1],uinf[n,2],beta,J);
            // Ts_obs[n,1] = t12[1,1];//gamma_rng(1/cvd^2,1/cvd^2/max([t12[1,1],0.5]));
            // Ts_obs[n,2] = t12[1,2];//gamma_rng(1/cvd^2,1/cvd^2/max([t12[1,2],0.5]));
            Ts_obs[n,1] = gamma_rng(1/cvd^2,1/cvd^2/t12[1,1]);
            Ts_obs[n,2] = gamma_rng(1/cvd^2,1/cvd^2/t12[1,2]);
        }
    // lambda = mean(calcTdyn(J,beta));
    // for(n in 1:Nd)
    //     {
    //        Tdyn_obs[Nd] = poisson_rng(lambda);
    //     } 
   
    }



