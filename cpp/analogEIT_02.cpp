/*
 * File: analogEIT_02.cpp
 *
 * This program integrates...
 *
 * It requires the following files to compile:
 *   nr3.h, odeint.h, stepper.h, stepperdopr853.h
 *
 * Use the following command line sequence to compile this file:
 *   g++ -lm analogEIT_01.cpp
 *
 *  Comments added by Juan D. Serna - October 12, 2017
 */

#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
//#include "rk4.h"
//#include "stepperdopr5.h"
#include "stepperdopr853.h"

#include <cmath>
#include <fstream>

// Initial conditions
Doub m1 = 2.56;
Doub m2 = 2.72;
Doub k1 = 400.0;
Doub k2 = 400.0;
Doub w1 = sqrt(k1/m1);
Doub w2 = sqrt(k1/m2);
Doub w = w1-0.0;
Doub K = 200.0;
Doub gm1 = 2.42;  //4.0E-2;
Doub gm2 = 0.54;  //1.0E-7;
Doub Omg1 = sqrt(K/m1);
Doub Omg2 = sqrt(K/m2);
Doub F = 1.0;
Doub C = 0.1;
Doub g = 9.8;

Doub ti = 0.0;
Doub tf = 500.0;
Doub x1i = 0.0;
Doub v1i = 0.0;
Doub x2i = 0.0;
Doub v2i = 0.0;

Int nsteps = 10000;  // Number of data points to print out

// Functor defining the RHS derivatives and the Jacobi matrix
struct rhs_EIT
{
    Doub dummy1;
    rhs_EIT(Doub dummy2) : dummy1(dummy2) {}

    // User-supplied routine that computes the derivatives of the RHS with
    // respect to x
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx)
    {
        Doub Coulomb1;
        Doub Coulomb2;

        if (y[1] > 0) { Coulomb1 = -C*m1*g; }
        else if (y[1] < 0) { Coulomb1 = C*m1*g; }
        else { Coulomb1 = 0; }

        if (y[5] > 0) { Coulomb2 = -C*m2*g; }
        else if (y[5] < 0) { Coulomb2 = C*m2*g; }
        else { Coulomb2 = 0; }

        dydx[0] = y[1];
        dydx[1] = -gm1*y[1] - w1*w1*y[0] + Omg1*Omg1*y[4] + (F/m1)*cos(w*x) + Coulomb1;
        dydx[2] = y[3];
        dydx[3] = -gm1*y[3] - w1*w1*y[2] + Omg1*Omg1*y[6] - (F/m1)*sin(w*x);
        dydx[4] = y[5];
        dydx[5] = -gm2*y[5] - w2*w2*y[4] + Omg2*Omg2*y[0] + Coulomb2;
        dydx[6] = y[7];
        dydx[7] = -gm2*y[7] - w2*w2*y[6] + Omg2*Omg2*y[2];

        //cout << fixed << setprecision(6);
        //cout << setw(12) << Coulomb1 << setw(12) << y[1]
             //<< setw(12) << Coulomb2 << setw(12) << y[5] << endl;
    }
};

// Main function
Int main(void)
{
    // Algorithm parameters
    const Int nvar = 8;        // Number of variables (equations to integrate)
    const Doub atol = 1.0E-3;  // Absolute tolerance (error)
    const Doub rtol = atol;    // Relative tolerance (rtol = atol works good)
    const Doub h1 = 0.01;      // Initial guess for a stepsize
    const Doub hmin = 0.0;     // Minimum stepsize (could be zero)

    // Some variable's declaration
    const Doub x1 = ti;
    const Doub x2 = tf;
    VecDoub ReP1(nsteps), ReP2(nsteps);
    VecDoub ystart(nvar);

    // Instantiation of initial conditions
    ystart[0] = x1i;
    ystart[1] = v1i;
    ystart[2] = 0.0;
    ystart[3] = 0.0;
    ystart[4] = x2i;
    ystart[5] = v2i;
    ystart[6] = 0.0;
    ystart[7] = 0.0;

    // Number of integrated points in the interval
    // (MAXSTP=50000 may be changed in odeint.h)
    Output out(nsteps);  // Dense output at 20 points plus x1
    rhs_EIT EITeqns(C);   // Declare EITeqns as an rhs_EIT object.

    // Integration routine
    Odeint < StepperDopr853 <rhs_EIT> >
    //Odeint < StepperDopr5 <rhs_EIT> >
    //Odeint < rk4 <rhs_EIT> >
    ode(ystart, x1, x2, atol, rtol, h1, hmin, out, EITeqns);
    ode.integrate();

    // Get the Re(Power) dissipated by m1 and m2
    for(Int i = 0; i < out.count; i++)
    {
        Doub tt = out.xsave[i];

        ReP1[i] = F*(out.ysave[1][i]*cos(w*tt) + out.ysave[3][i]*sin(w*tt));
        ReP2[i] = F*(out.ysave[5][i]*cos(w*tt) + out.ysave[7][i]*sin(w*tt));
    }

    // Print output
    cout << fixed << setprecision(6);
    for(Int i = 0; i < out.count; i++)
    {
        cout << setw(12) << out.xsave[i]
             << setw(12) << out.ysave[0][i]
             << setw(12) << out.ysave[1][i]
             << setw(12) << out.ysave[2][i]
             << setw(12) << out.ysave[3][i]
             << setw(12) << out.ysave[4][i]
             << setw(12) << out.ysave[5][i]
             << setw(12) << out.ysave[6][i]
             << setw(12) << out.ysave[7][i]
             << setw(12) << ReP1[i]
             << setw(12) << ReP2[i]
             << endl;
    }

    //cout << endl;
    //cout << ttime.size() << "\t"
    //<< vm1.size() << "\t"
    //<< Pm1.size() << "\t"
    //<< out.ysave.nrows() << "\t"
    //<< endl;

    return 0;
}
