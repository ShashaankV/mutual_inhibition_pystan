// fits static rivalry and intermittent (dynamic S)
// assume shared parameters except gamma
// make sure tau are consistent, tau_a in Wilson for static is unit magnitude
// dynamic S has tau_u and tau_a, so tau_u is set to tau_a/1000
// dt is tau_u/10

functions
{
#include ./functions/wilson_fun.stan
#include ./functions/dynamicS_fun.stan
}

data 
    {
        int<lower = 0> Ns; // number of static rivalry measurements
        real Ts_obs[Ns, 2]; // measured T static 
        real c[Ns,2]; //contrast (stimulus) condition
        real x_r[1]; // junk for solver
        int x_i[1]; // junk
        int max_steps;
        real rel_tol;
        real function_tol;

        int<lower = 0> Nd; // number of intermittent rivalry measurements
                             // observation, report noise;
        real Tdyn_obs[Nd]; //noisy report during intermittent rivalry 
    }

transformed data
    {
        real cvd;
        cvd = 0.6;
    }

parameters 
    {
        real<lower=0> J;
        real<lower=0> beta;
    }



transformed parameters
    { 
        //static rivalry variables
        real Ts_est[Ns, 2]; // need to reformat this so only calc each c-set once
        real uinf[Ns,2];
        real t12[1,2];
        real<lower=1e-10, upper=1> u1logarg[Ns]; // not allowed to sample out of wilson rivalry
        real<lower=1e-10, upper=1> u2logarg[Ns];
        //intermittent rivalry variables 
        real Tdyn_est; //estimate of the mean for single fixed condition

        //calc static
        for(n in 1:Ns)
        {
            vector[2] uinf_input_1 = [c[n,1],J]';
            vector[2] uinf_input_2 = [c[n,2],J]';
            vector[1] guess = [0.01]';
            vector[1] sol_1 = algebra_solver(f, guess, uinf_input_1, x_r, x_i,rel_tol,function_tol,max_steps);
            vector[1] sol_2 = algebra_solver(f, guess, uinf_input_2, x_r, x_i,rel_tol,function_tol,max_steps);
            uinf[n,1] = sol_1[1];
            uinf[n,2] = sol_2[1];

            u1logarg[n] = (J*c[n,1]^0.1 - beta*uinf[n,2])/(uinf[n,2]);
            u2logarg[n] = (J*c[n,2]^0.1 - beta*uinf[n,1])/(uinf[n,1]);

            t12 = T12fun(c[n,1],c[n,2],uinf[n,1],uinf[n,2],beta,J);
            Ts_est[n,1] = t12[1,1];
            Ts_est[n,2] = t12[1,2];
        } 

        //calc intermittent
        Tdyn_est = calcTdyn(J,beta);
         
    }

model
    {
        beta ~ cauchy(1.,1.); 
        J ~ cauchy(1.,1.); 

        for(n in 1:Ns)
        {
            target += gamma_lpdf(Ts_obs[n,1]|1/cvd^2,1/cvd^2/max([Ts_est[n,1],0.5]));
            target += gamma_lpdf(Ts_obs[n,2]|1/cvd^2,1/cvd^2/max([Ts_est[n,2],0.5]));

            u1logarg[n] ~ uniform(1e-10,1);
            u2logarg[n] ~ uniform(1e-10,1);
        }

        for(n in 1:Nd)
        {
            target += poisson_lpmf(Tdyn_obs[n]|Tdyn_est); // maybe gamma
        }
    }
generated quantities //nothing here for backward other than ll
    {    
    }



