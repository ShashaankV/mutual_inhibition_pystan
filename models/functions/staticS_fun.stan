real[] calcTs(real dt,int iterations,real[,] S)
        {
            real gamma;
            real beta;
            real J;
            real tau_u;
            real tau_a;
            real a[2];
            real u[2,iterations];
            real x;
            
            int st; //current dominance state
            int st_next; //next dominance state
            int count; //counter for current epoch duration
            real T1; //epoch duration 
            real T2; //epoch duration 
            int epochs1; //counter for percept epochs
            int epochs2; //counter for percept epochs
            real T[2];
            
            T1 = 0;
            T2 = 0;
                        
            gamma = 1.8;
            beta = 2;
            J = 1;
            tau_u = 1;
            tau_a = 1000;
            
            a[1] = 0;
            a[2] = 0;
            
            u[1,1] = 0.01;
            u[2,1] = 0.1;
            
            epochs1 = 0;
            epochs2 = 0;
            
            if(u[1,1]>u[2,1]){st = 0;}
            else{st = 1;}
            
            count = 0;
            
            for(it in 2:iterations)
            {
               count = count +1;
               x = S[1,it] - beta*u[2,it-1] - gamma*a[1];
               u[1,it] = u[1,it-1] + dt/tau_u*(-u[1,it-1] + g(x));
               a[1] = a[1] + dt/tau_a*(-a[1] + u[1,it]);
               x = S[2,it] - beta*u[1,it-1] - gamma*a[2];
               u[2,it] = u[2,it-1] + dt/tau_u*(-u[2,it-1] + g(x));
               a[2] = a[2] + dt/tau_a*(-a[2] + u[2,it]);
               
               
              if(u[1,it]>u[2,it]){st_next = 0;}else{st_next = 1;}
                
              if(st!=st_next){
              if(st==0){T1 = (T1 + count*dt);epochs1 = epochs1 + 1;}
                   else{T2 = (T2 + count*dt);epochs2 = epochs2 + 1;}
              st = st_next;
              count = 0;
              }
              
             } 
            T = {T1/epochs1,T2/epochs2};
            return T;
        } 