// gamma, alpha, tauu, taua, etag, etaS, c0 are hardcoded and fixed
// small constant drive is c0

//uinffun and f are used for algebra solver to get uinf

real uinffun(real z, real c, real J)
    {
        real uinf;
        real x;
        real gamma;

        gamma = 2;
        x = J*c^0.1 + 1e-3 - gamma*z; //implicit alpha = 0 and gamma = -2
        uinf = ((x>0)*x+1e-6)^0.5;
        return uinf;
    }

vector f(vector uinf, vector theta, real[] x_r, int[] x_i)
        {
            return [uinf[1] - uinffun(uinf[1], theta[1], theta[2])]';
        }

real[,] T12fun(real c1, real c2, real u1, real u2, real beta, real J)
    {
        //u1 and u2 calculated in transformed parameters
        real T[1,2];
        real tau_a;
        real x1;
        real x2;
        real gamma;

        gamma = 2;

        tau_a = 1;

        x2 = J*c1^0.1 + 1e-3 - beta*u2;
        x1 = J*c2^0.1 + 1e-3 - beta*u1;
        //wta condition (if beta++, then possible for x1 and x2 <0, may need to set one >0)
        if (x1 <= 0 || x2 <=0)
            {
                T[1,1] = (x1 > 0)*60 + 1e-10; // if x1 can be wta then set T1 to block duration, 60 seconds
                T[1,2] = (x2 > 0)*60 + 1e-10;
            }
        //normalization condition
        else if (u1<=0 || u2<=0 || u1>1000 || u2>1000)
        {
            T[1,1] = 1e-10;
            T[1,2] = 1e-10;
        }
        else if (x1/(gamma*u2)>1 || x2/(gamma*u1)>1)
            {
                T[1,1] = 1e-10;
                T[1,2] = 1e-10;
            }
        else
            {
                T[1,1] = -tau_a*log(x1/(gamma*u2)); //implicit gamma = -1
                T[1,2] = -tau_a*log(x2/(gamma*u1));
            }
        return T;
        }

