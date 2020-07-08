// synthesizes static rivalry using Wilson approximation

functions
{
#include ../models/functions/dynamicS_fun.stan
}

data 
    {
        real J; 
        real beta;
        int iterations;
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
        real lambda;
        int Tdyn_obs[Nd]; //synthetic noisy report during intermittent rivalry 
        real u[2,iterations];
        
        lambda = mean(calcTdyn(J,beta,iterations));
        for(n in 1:Nd)
            {
            Tdyn_obs[n] = poisson_rng(lambda);
            }  
        u = calcUdyn(J,beta,iterations);
    }



