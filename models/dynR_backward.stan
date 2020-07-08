// synthesizes static rivalry using Wilson approximation

functions
    {
    real [,] calcU(
            real t,             // time
            real[] u,           // system state variables: u1,u2,a1,a2
            real[] theta,       // parameters
            real[] x_r,
            int[] x_i
            )
              
              real J = theta[1];
              real beta = theta[2];
              
              real c0 = 1e-4;
              real gamma = 1e-3;
              real T = 500;
              real Ton = 50;
              real tau_a = 1000;
            
              real u1 = u[1];
              real u2 = u[2];
              real a1 = u[3];
              real a2 = u[4];
              real x;
              
              x = J*(t%Ton<T) + c0 - beta*u2 - gamma*a1;
              du_dt[1] = -u1 + ((x<0)*x)^0.5; // u1
              
              x = J*(t%Ton<T) + c0 - beta*u1 - gamma*a2;
              du_dt[2] = -u2 + ((x<0)*x)^0.5; // u2
              
              du_dt[3] = (-a1 + u1)/tau_a; // a1
              du_dt[4] = (-a2 + u2)/tau_a; // a2
              
              return du_dt;    
        }

data 
    {
        int<lower = 0> Nd; // number of intermittent rivalry measurements
                             // observation, report noise;
        int u_obs[Nd]; //synthetic noisy report during intermittent rivalry
    }


parameters
    {
        real<lower=0, upper=10> J; 
        real<lower=0, upper=10> beta;
    }

transformed parameters
    {
        real u[2,iterations];
        real<lower=1> lambda;
        lambda = mean(calcTdyn(J,beta));
    }

model //
    {   

        beta ~ cauchy(1.,1.); 
        J ~ cauchy(1.,1.); 

        for(n in 1:Nd)
            {
            target += poisson_lpmf(Tdyn_obs[n]|lambda);
            }  
    }



