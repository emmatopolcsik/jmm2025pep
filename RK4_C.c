// Runge Kutta 4th Order
#include <math.h> /*Include library */
#include "mex.h" /*Always include this */

// EOM with STM
void eomSTM_C(double t,double *x,double mu,double *dx)
{
double xx = x[0];
double yy = x[1];
double zz = x[2];
double vx = x[3];
double vy = x[4];
double vz = x[5];

double yy2pzz2 = pow(yy,2) + pow(zz,2);
double xpmu = xx+mu;
double xpmu2 = pow(xpmu,2);
double xpmum1 = xpmu-1;
double xpmum12 = pow(xpmum1,2);
double d2 = xpmu2 + yy2pzz2;
double r2 = xpmum12 + yy2pzz2;
double muoverr3 = mu/pow(r2,1.5);
double onemmuoverd3 = (1-mu)/pow(d2,1.5);
double monemmuoverd3mmuoverr3 = -(onemmuoverd3 + muoverr3);
double monemmuoverd3mmuoverr3p1 = monemmuoverd3mmuoverr3 + 1;
double threemuoverr5 = 3*mu/pow(r2,2.5);
double threeonemmuoverd5 = 3*(1-mu)/pow(d2,2.5);
double threeonemmuoverd5p3muoverr5 = threeonemmuoverd5+threemuoverr5;

double Uxx, Uyy, Uzz, Uxy, Uxz, Uyz;
Uxx = monemmuoverd3mmuoverr3p1 + threeonemmuoverd5*xpmu2 + threemuoverr5*xpmum12;
Uyy = monemmuoverd3mmuoverr3p1 + pow(yy,2)*threeonemmuoverd5p3muoverr5;
Uzz = monemmuoverd3mmuoverr3   + pow(zz,2)*threeonemmuoverd5p3muoverr5;
Uxy = threeonemmuoverd5*xpmu*yy + threemuoverr5*xpmum1*yy;
Uxz = threeonemmuoverd5*xpmu*zz + threemuoverr5*xpmum1*zz;
Uyz = (yy*zz)*threeonemmuoverd5p3muoverr5;

double A[6][6] =
{
    {0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 1},
    {Uxx, Uxy, Uxz,  0, 2, 0},
    {Uxy, Uyy, Uyz, -2, 0, 0},
    {Uxz, Uyz, Uzz,  0, 0, 0},
};

int i, j, k;
double phi1i, phi2i, phi3i, phi4i, phi5i, phi6i;
int idx0, idx1, idx2, idx3, idx4, idx5;
for(i=1; i<7; i++){ 
    idx0 = 6*i;
    idx1 = idx0+1;
    idx2 = idx0+2;
    idx3 = idx0+3;
    idx4 = idx0+4;
    idx5 = idx0+5;
    phi1i = x[idx0];
    phi2i = x[idx1];
    phi3i = x[idx2];
    phi4i = x[idx3];
    phi5i = x[idx4];
    phi6i = x[idx5];
    dx[idx0] = phi4i;
	dx[idx1] = phi5i;
	dx[idx2] = phi6i;
    j=3;
    dx[idx3] = A[j][0]*phi1i+A[j][1]*phi2i+A[j][2]*phi3i+A[j][4]*phi5i;
    j=4;
    dx[idx4] = A[j][0]*phi1i+A[j][1]*phi2i+A[j][2]*phi3i+A[j][3]*phi4i;
    j=5;
    dx[idx5] = A[j][0]*phi1i+A[j][1]*phi2i+A[j][2]*phi3i;
}

dx[0] = vx;
dx[1] = vy;
dx[2] = vz;
dx[3] = -onemmuoverd3*xpmu - muoverr3*xpmum1 + xx + 2*vy;
dx[4] = yy*monemmuoverd3mmuoverr3p1 - 2*vx;
dx[5] = zz*monemmuoverd3mmuoverr3;
}


// Runge-Kutta 4th Order Integrator
// Inputs: dt (time elapsed) 
//         x0[42] (initial state)
//         mu (parameter)
//         steps (number of steps in propagation)
//         xn[42] (output state after dt)
void RK4(double dt, double *x0, double mu, double steps, double *xn){
    double h = dt/steps;

    // Define xn = x0
    for (int i = 0; i < 42; i++){
        xn[i] = x0[i];
    }
    
    double a[42];
    double b[42], opB[42];
    double c[42], opC[42];
    double d[42], opD[42];
    double val[42], xn_1[42];

    // Enter iterative loop
    for (int j = 0; j < steps; j++){
        // a = f(y0,t0) = f(xn)
        eomSTM_C(0,xn,mu,a);

        // b = f(y0 + a*h/2) = f(xn + a*h/2)
        for (int i = 0; i < 42; i++){
            opB[i] = xn[i] + (h/2)*a[i];
        }
        eomSTM_C(0,opB,mu,b);

        // c = f(y0 + b*h/2) = f(xn + b*h/2)
        for (int i = 0; i < 42; i++){
            opC[i] = xn[i] + (h/2)*b[i];
        }
        eomSTM_C(0,opC,mu,c);

        // d = f(y0 + c*h) = f(xn + c*h)
        for (int i = 0; i < 42; i++){
            opD[i] = xn[i] + h*c[i];
        }
        eomSTM_C(0,opD,mu,d);

        // define next time step = xn + h/6(a + 2b + 3c + d)
        for (int i = 0; i < 42; i++){
            val[i] = h/6*(a[i] + 2*b[i] + 2*c[i] + d[i]);
            xn_1[i] = xn[i] + val[i];
        }

        // prepare for next loop
        //mexPrintf("[%i] Step: ",j);
        for (int i = 0; i < 42; i++){
            xn[i] = xn_1[i];
            //if (i<6){
            //mexPrintf("%f ", xn[i]);}
            //}
            //mexPrintf("\n");
        }
    }
    
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
// Check the number of input and output arguments
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidNumInputs", "Invalid number of input arguments. Expected 4.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidNumOutputs", "Invalid number of output arguments. Expected 1.");
    }
    // Initialize inputs/outputs
    double mu, dt, steps;
    double *x0, *xn;
    // Get the value of the scalar inputs
    dt           = mxGetScalar(prhs[0]);
    mu           = mxGetScalar(prhs[2]);
    steps        = mxGetScalar(prhs[3]);

    // Create a pointer to the real data in the input vectors/matrices
    x0        = mxGetPr(prhs[1]);

    // Creating output matrices
    plhs[0] = mxCreateDoubleMatrix(42, 1, mxREAL);

    // Get a pointer to the real data in the output matrix
    xn = mxGetPr(plhs[0]);

    // Call the computational function
    RK4(dt, x0, mu, steps, xn);

}