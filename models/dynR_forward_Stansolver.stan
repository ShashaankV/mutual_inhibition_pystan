// synthesizes static rivalry using Wilson approximation

functions
    {
real [] calcU(
            real t,             // time
            real[] u,           // system state variables: u1,u2,a1,a2
            real[] theta,       // parameters
            real[] x_r,
            int[] x_i
            )
         { 
              real du_dt[4];
              real J = theta[1];
              real beta = theta[2];
              
              real c0 = 1e-4;
              real gamma = 1e-3;
              real T = 500;
              int Ton = 50;
              real tau_a = 1000;
            
              real u1 = u[1];
              real u2 = u[2];
              real a1 = u[3];
              real a2 = u[4];
              real x;
              
              real mod;
              
              mod = t - T*floor(t/T);
              
              x = J*(mod<Ton) + c0 - beta*u2 - gamma*a1;
              du_dt[1] = -u1 + ((x>0)*x)^0.5; // u1
              
              x = J*(mod<Ton) + c0 - beta*u1 - gamma*a2;
              du_dt[2] = -u2 + ((x>0)*x)^0.5; // u2
              
              du_dt[3] = (-a1 + u1)/tau_a; // a1
              du_dt[4] = (-a2 + u2)/tau_a; // a2
              
              return du_dt; 
        } 
        
    real calcT(real [] x1, real[] x2, int Nts)
        {
            int k;
            real counter;
            real dx;
            real T = 0;
            real s;
            
            dx = (x1[1] - x2[1])>0;
            s = dx;
            counter = 1;
            k = 0;
            for(i in 2:Nts){
                counter += 1;
                dx = (x1[i] - x2[i])>0;
                if(s!=dx){k += 1; T = T + counter; s = dx; counter = 0;}
                }
            if(k>0){T=T/k;}
            else{T=1000.0;}
            return T;
         }

    }

data 
    {
        int<lower = 0> Nts; // number of milliseconds
        real t[Nts]; //
        real t0;
        real J;
        real beta;
        real x_r[1];
        int x_i[1];
        int Nd;
    }

transformed data
    {
        real cvd;
        cvd = 0.6;
    }


generated quantities //
    {   
        real theta[2] = {J,beta};
        real u[Nts, 4];   // solution from the ODE solver
        real u_init[4];
        real T;
        real Tdyn_obs[Nd];

        u_init[1] = 0.1;
        u_init[2] = 0.01;
        u_init[3] = 0;
        u_init[4] = 0;

        u = integrate_ode_rk45(calcU, u_init, t0, t, theta, x_r, x_i);//,1e-2,1e-2,2000);  
        T = calcT(u[:,3],u[:,4],Nts);
        
        for(n in 1:Nd){
            Tdyn_obs[n] = gamma_rng(1/cvd^2,1/cvd^2/T);
         }
    }



