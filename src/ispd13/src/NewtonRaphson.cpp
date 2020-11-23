/// \file NewtonRaphson.cpp
/// This is a demonstration of how to use NewtonRaphson.h<br />
/// to compile:  g++ -I. -lm NewtonRaphson.cpp<br />
/// description: this test the algorithm in 3 situations:<br />
///    1) calculating df at bouth points, to find the average df, for each iteration, with no boundaries <br />
///    2) calculating d2f to get second order approximation at each iteration, with no boundaries<br />
///    3) calculating d2f to get second order approximation at each iteration, and a boundary solution<br />

/***************************************************************************
 *   Copyright (C) 2009 by Clark Sims                                      *
 *   http://AcumenSoftwareInc.com/WhoWeAre/Clark_Sims.html                 *
 *   ClarkSims@AcumenSoftwareInc.com                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "NewtonRaphson.h"
#include <math.h>

/// \brief This functor is provided as a demonstration of how to create a functor for the template class NetwonRaphson0. It encapsulates the function: exp(x) - offset.  The function variable is x. Offset is a hidden variable. This class is defined in the demonstration file NewtonRaphson.cpp.
class Newton_functor {
public:
  /// This is the hidden variable for class Newton_functor, which makes this class a functor, as opposed to a group of function pointers.
    double Rd;
    double R;
    double C1;
    double C2;
    double dt;

  /// the creator for Newton_functor which assigns the hidden variable offset
  Newton_functor( double Rd, double R, double C1, double C2, double dt  ) : Rd(Rd), R(R), C1(C1), C2(C2), dt(dt) { }

  /// the function to be solved
  double f( double Ceff) {
	  //2.282849659*10^(-22)-6.00*10^(-10)*Ceff+9*Ceff^2*(1-exp(-6.666666667*10^(-11)/Ceff))
    //return 2.282849659e-22 - (6.00e-10 * x) + 9*x*x*(1-exp(-6.666666667e-11/x));
      
      double z = (C2+C1)/(C2*R*C1);
      double a = Rd*R*C1*C2;
      double b = (C1 + C2)*Rd + R*C1;
      double inRoot = b*b - 4.0*a;
      double p1 = (b + sqrt(inRoot))/(2.0*a);
      double p2 = (b - sqrt(inRoot))/(2.0*a);
      
      double A = z/(p1*p2);
      double B = (z-p1)/(p1*(p1-p2));
      double D = (z-p2)/(p2*(p2-p1));
      
      double RdCeff = Rd*Ceff;
      double p = dt/(RdCeff);
      double Rddt = Rd*dt;
      
      return A*dt + (B/p1)*(1-exp(-p1*dt)) + (D/p2)*(1-exp(-p2*dt)) - RdCeff*dt + RdCeff*RdCeff*(1-exp(-p));
  }

  /// the derivative of the function to be solved
  double df( double Ceff) {
      double RdCeff = Rd*Ceff;
      double p = dt/(RdCeff);
      double Rddt = Rd*dt;
      return -Rddt + 2.0*Rd*RdCeff*(1-exp(-p)) - Rddt*exp(-p);
  }
};

int NRmain() {
	double Rd = 3000.0, R = 400.0, C1 = 0.500e-12, C2 = 0.200e-12, dt = 0.200e-9;
    Newton_functor functor( Rd, R, C1, C2, dt);
    // {
    NewtonRaphsonSolve0< Newton_functor, double >
    newton(
           1.2e-13,               // X0
           1e-30,                 // Epsilon
           true,                 // Twice_df
           100,                  // max_iter
           functor,              // F
           &Newton_functor::f,   // f
           &Newton_functor::df); // df
    
    //   cout << "running iteration with f and df defined, over all real numbers" << endl;
    newton.set_check_boundary( true);
    newton.set_max_x( C1 + C2 );
    newton.set_min_x( C2 );
    cout << "result: " << newton.do_iteration( NULL/*&cout*/);
    cout << endl;

/*
    cout << "running iteration with f, df and d2f defined, over all real numbers" << endl;
    newton.set_d2f( &Newton_functor::df);
    newton.do_iteration( &cout);

    cout << "running iteration with f, df and d2f defined, over the interval [1.6, infinity)" << endl;
    newton.set_check_boundary( true);
    newton.set_max_x( 1.6);
    newton.do_iteration( &cout);
  }

  cout << "running a numerically equivalent example, where functor._offset = 0, and newton._fsolve = 5.0"
       << endl;

  functor._offset = 0;
  {
    NewtonRaphsonSolve0< Newton_functor, double >
      newton(
	     1.2,                  // X0
	     1e-8,                 // Epsilon
	     true,                 // Twice_df
	     100,                  // max_iter
	     functor,              // F
	     &Newton_functor::f,   // f
	     &Newton_functor::df,  // df
	     0,                    // d2f
	     5.0);                 // Fsolve

    cout << "running iteration with f and df defined, over all real numbers" << endl;
    newton.do_iteration( &cout);

    cout << "running iteration with f, df and d2f defined, over all real numbers" << endl;
    newton.set_d2f( &Newton_functor::df);
    newton.do_iteration( &cout);

    cout << "running iteration with f, df and d2f defined, over the interval [1.6, infinity)" << endl;
    newton.set_check_boundary( true);
    newton.set_max_x( 1.6);
    newton.do_iteration( &cout);
  }
*/
  return 0;
}

/*! \mainpage NewtonRaphson

<p>
This library contains a template function NewtonRaphsonSolve0 which is an implentation of the Newton Raphson algorithm. The template uses two type, functor, and real. "real" can be any implentation of one dimensional floating point arithmetic. "functor" is a class, wich must have a member function pointer, "_f", which will be the functions which is solveved. There are several other arguments, which can be used to tweek the algorithm. These arguments are described in the documentation for the constructor of the template class NewtonRaphsonSolve0::NewtonRaphsonSolve0. The actual iteration is performed by the member function NewtonRaphsonSolve0::do_iteration.
</p>

<p>
One can view or download the original source code for NewtonRaphson.h <a href="NewtonRaphson.h"> here </a>. This header file contains the code which defines the template class NewtonRaphsonSolve0. <br/>
One can view or download the original source code for NewtonRaphson.cpp <a href="NewtonRaphson.cpp"> here </a>. This source file contains three examples which call call objects of class NewtonRaphsonSolve0.
</p>

<p>
This computer code is being released under the GNU general public license:
</p>

<pre>
Copyright (C) 2009 by Clark Sims
http://AcumenSoftwareInc.com/WhoWeAre/Clark_Sims.html
ClarkSims@AcumenSoftwareInc.com

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the
Free Software Foundation, Inc.,
59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
</pre>

*/
