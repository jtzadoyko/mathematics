/*
 * File: Lorenz.cpp
 *
 * This program integrates a stiff set of differential equations (Lorenz
 * Equations) using Rosenbrock and Semi-Implicit Extrapolation methods.
 *
 * It requires the following files to work:
 *
 * nr3.h, odeint.h, ludcmp.h, stepper.h, stepperross.h, <steppersie.h>
 *
 * Use the following command line sequence to compile this file:
 *    g++ -lm Lorenz.cpp
 *
 *  Comments added by Juan D. Serna - October 8, 2017
 */

#include "nr3.h"
#include "odeint.h"
#include "ludcmp.h"
#include "stepper.h"
#include "stepperross.h"  // Use this for Rosenbrock methods
//#include "steppersie.h"   // Use this for Semi-Implicit Extrapolation

#include <cmath>
#include <fstream>

// Lorenz parameters
double BETA = 8.0/3.0;
double RHO = 28.0;
double SIGMA = 10.0;

// Functor defining the RHS derivatives and the Jacobi matrix
struct rhs_Lorenz
{
    Doub dummy;
    rhs_Lorenz(Doub ddummy) : dummy(ddummy) {}

    // User-supplied routine that computes the derivatives of the RHS with
    // respect to x
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx)
    {
        // Lorenz equations
        dydx[0] = SIGMA*(y[1] - y[0]);
        dydx[1] = y[0]*(RHO - y[2]) - y[1];
        dydx[2] = y[0]*y[1] - BETA*y[2];
    }

    // User-supplied routine that computes the Jacobi matrix of derivatives of
    // the RHS with respect to the components of y
    void jacobian(const Doub x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy)
    {
        Int n = y.size();

        // explicit dependence on x
        for (Int i = 0; i<n; i++)
          dfdx[i] = 0.0;

        // Jacobian matrix elements
        dfdy[0][0] = -SIGMA;
        dfdy[0][1] = SIGMA;
        dfdy[0][2] = 0.0;
        dfdy[1][0] = RHO - y[2];
        dfdy[1][1] = -1.0;
        dfdy[1][2] = -y[0];
        dfdy[2][0] = y[1];
        dfdy[2][1] = y[0];
        dfdy[2][2] = -BETA;
    }
};

// Main function
Int main(void)
{
    // Algorithm's initial conditions
    const Int nvar = 3;         // Number of variables (equations to integrate)
    const Doub atol = 1.0E-10;  // Absolute tolerance (error)
    const Doub rtol = atol;     // Relative tolerance (rtol = atol works good)
    const Doub h1 = 1.0E-3;     // Initial guess for a stepsize
    const Doub hmin = 0.0;      // Minimum stepsize (could be zero)

    // Initial conditions of the problem to integrate
    VecDoub ystart(nvar);
    const Doub x1 = 0.0;   // Lower bound of integration
    const Doub x2 = 50.0;  // Upper bound of integration
    ystart[0] = 1.0;
    ystart[1] = 1.0;
    ystart[2] = 1.0;

    // Number of integrated points in the interval
    // (MAXSTP=50000 may be changed in odeint.h)
    Output out(10000);     // Dense output at 20 points plus x1
    rhs_Lorenz d(1.0E-3);  // Declare d as a rhs_Lorenz object.

    // Integration routine
    //Odeint < StepperSie <rhs_Lorenz> >
        //ode(ystart, x1, x2, atol, rtol, h1, hmin, out, d);
    Odeint < StepperRoss <rhs_Lorenz> >
        ode(ystart, x1, x2, atol, rtol, h1, hmin, out, d);
    ode.integrate();

    // Print output
    cout << fixed << setprecision(6);
    for(Int i = 0; i < out.count; i++)
    {
        cout << setw(12) << out.xsave[i]
             << setw(12) << out.ysave[0][i]
             << setw(12) << out.ysave[1][i]
             << setw(12) << out.ysave[2][i]
             << endl;
    }

    return 0;
}
