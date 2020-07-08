// dynamic rivalry dominace duration (measured in pings), 
//expects S with pings every Tp with duty cycle Tpdur/Tp

// fix Son base, Soff, tau's, and other time stuff
// J scales Son base 

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
            k = -1;
            for(i in 2:Nts){
                counter += 1;
                dx = (x1[i] - x2[i])>0;
                if(s!=dx){k += 1; T = T + counter; s = dx; counter = 0;}
                }
            return T/k;
         }


// real g(real x)
//     {
//         real y;
//         y = ((x>0)*x)^0.5; //gain and threshold are fixed
//         return y;
//     }
      
// // need to fix so not saving u-series in loop
//  real[] calcTdyn(real J,real beta)
//         {

//             int Tp;
//             real Tpdur;
//             real dt;

//             real Son;
//             real Soff;
//             real gamma;
//             // real beta;
//             // real J;
//             real tau_u;
//             real tau_a;
//             real a[2];
//             real u[2,iterations];
//             real x;

//             int st; //current dominance state
//             int st_next; //next dominance state
//             int count; //counter for current epoch duration (unit = pings)
//             real T1; //epoch duration (unit = pings)
//             real T2; //epoch duration 
//             int epochs1; //counter for percept epochs (for averaging), 
//                          //deterministic but not keeping track of multiple epochs in output 
//             int epochs2; //counter for percept epochs
//             real T[2];
            
//             Son = 200;
//             Soff = 1e-4;
            
//             Tp = 5000; //ping every Tp iterations
//             Tpdur = 1; 
//             dt = 1/10.0;
                                  
//             gamma = 1e-3;
//             // beta = 2;
//             // J = 1;
//             tau_u = 1;
//             tau_a = 1000;
            
//             a[1] = 0;
//             a[2] = 0;
            
//             u[1,1] = 0.01;
//             u[2,1] = 0.1;
            
//             T1 = 0;
//             T2 = 0;
            
//             epochs1 = 0;
//             epochs2 = 0;
            
//             if(u[1,1]>u[2,1]){st = 0;}
//             else{st = 1;}
            
//             count = 0;

//             for(it in 2:iterations)
//             {
//                count = count + (it%Tp==0);
               
//                x = J*(it%Tp<Tpdur)*Son + Soff - beta*u[2,it-1] - gamma*a[1];
               
//                u[1,it] = u[1,it-1] + dt/tau_u*(-u[1,it-1] + g(x));
//                a[1] = a[1] + dt/tau_a*(-a[1] + u[1,it]);
               
//                x = J*(it%Tp<Tpdur)*Son + Soff - beta*u[1,it-1] - gamma*a[2];
               
//                u[2,it] = u[2,it-1] + dt/tau_u*(-u[2,it-1] + g(x));
//                a[2] = a[2] + dt/tau_a*(-a[2] + u[2,it]);
               
//                //check state near end of ping
//                if(it%Tp==(Tpdur-1))
//                {
//                   if(u[1,it]>u[2,it]){st_next = 0;}else{st_next = 1;}

//                   if(st!=st_next){
//                   if(st==0){T1 = (T1 + count);epochs1 = epochs1 + 1;}
//                        else{T2 = (T2 + count);epochs2 = epochs2 + 1;}
//                   st = st_next;
//                   count = 0;
//                    }
//                  }
//                  }
//             if ((epochs1*epochs2) > 0){T = {T1/epochs1,T2/epochs2};}
//             else {T = {1,1};}
//             return T;
        
//         } 



// real[,] calcUdyn(real J,real beta)
//         {

//             int Tp;
//             real Tpdur;
//             real dt;
//             int iterations;

//             real Son;
//             real Soff;
//             real gamma;
//             // real beta;
//             // real J;
//             real tau_u;
//             real tau_a;
//             real a[2];
//             real u[2,iterations];
//             real x;
            
//             Son = 200;
//             Soff = 1e-4;
            
//             Tp = 5000; //ping every Tp iterations
//             Tpdur = 1; 
//             dt = 1/10.0;
                                  
//             gamma = 1e-3;
//             // beta = 2;
//             // J = 1;
//             tau_u = 1;
//             tau_a = 1000;
            
//             a[1] = 0;
//             a[2] = 0;
            
//             u[1,1] = 0.01;
//             u[2,1] = 0.1;
            


//             for(it in 2:iterations)
//             {
  
//                x = J*(it%Tp<Tpdur)*Son + Soff - beta*u[2,it-1] - gamma*a[1];
               
//                u[1,it] = u[1,it-1] + dt/tau_u*(-u[1,it-1] + g(x));
//                a[1] = a[1] + dt/tau_a*(-a[1] + u[1,it]);
               
//                x = J*(it%Tp<Tpdur)*Son + Soff - beta*u[1,it-1] - gamma*a[2];
               
//                u[2,it] = u[2,it-1] + dt/tau_u*(-u[2,it-1] + g(x));
//                a[2] = a[2] + dt/tau_a*(-a[2] + u[2,it]);
               

//             return u;
        
//         } 
        
//         }

