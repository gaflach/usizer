/// \file NewtonRaphson.h
///
/// \brief This is an encapsulation of the one dimensional NewtonRaphson iteration algorithm. It allows for calculating df twice at each iteration, or d2f, to increase the spead of convergance. It also can be used in situations with boundaries. The orignal source cand be viewed or downloaded <a href="NewtonRaphson.h">here</a>.

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
#ifndef NewtonRaphson_hpp
#define NewtonRaphson_hpp

#include <iostream>
#include <vector>

using namespace std;

/// \brief This template class is an encapsulation of the one dimensional NewtonRaphson iteration algorithm. It allows for calculating df twice at each iteration, or d2f, to increase the spead of convergance. It also can be used in situations with boundaries. The orignal source cand be viewed or downloaded <a href="NewtonRaphson.h">here.</a><br/>
///
///
/// An object of class NewtonRaphsonSolve0 is created with a call to the constructor, NewtonRaphsonSolve0::NewtonRaphsonSolve0, which takes arguments for the starting point, X0, the convergance stopping point, Epsilon, and a functor f, and a functor df, which is the derivative of f. <br/>
/// The function do_iteration, attempts to find the value of x, _x_converge, which solves the equation, (_F.*_f(_x_converge) - _f_solve) < _epsilon. <br/>
/// The second derivative d2f can optionally be supplied. If check_boundary is set to true, then boundary solutions can be solved also. The member function do_iteration, starts the iterative process. <br/>
///The file NewtonRaphson.cpp contains 3 examples of how this template class can be used. The source code for NewtonRaphson.cpp can be downloaded <a href="NewtonRaphson.cpp">here.</a><br/>
template<class functor, class real>
class NewtonRaphsonSolve0 {
 protected:
  // key attributes

  /// \brief The maximum number of iterations of the routine. See the documentation for do_iteration,
  /// to see how this value affects the algorithm.
  int _max_iter;

  /// \brief This is the starting point of the iteration. If _check_boundary is true, then
  /// _x0 must be between _min_x and _max_x.
  real _x0;

  /// \brief the convergence parameter. The calculations stop when F(x) < epsilon.
  real _epsilon;

  /// \brief the value of f to be solved. This is usually 0 in the academic literature.
  real _fsolve;

  /// \brief if true then the value of the derivative is calculated twice, to calculate the average value of df, over the interval dx. This is useful, if df is easier to calculate than f. For example dprice/dsigam, can be aproximated ver
  bool _twice_df;

  /// \brief this is the function whose 0 is being solved
  real (functor::*_f)(real);

  /// \brief this is the derivative
  real (functor::*_df)(real);

  /// \brief this is the second derivative. If set to null, then it isn't used
  real (functor::*_d2f)(real);

  /// \brief this is the object used to calculate the functions
  functor& _F;

  // functional attributes
  /// \brief if true, then the search for the zero, has boundaries
  bool _check_boundary;

  /// \brief this is the minimum value of the acceptible range
  real _min_x;

  /// \brief this is the maximum value of the acceptible range
  real _max_x;


  // derived attributes
  /// \brief the value of x such that f(x) = 0 +/- epsilon
  mutable real _x_converge;
  /// \brief the last value of f
  mutable real _f_converge;
  /// \brief the sequence of x's
  mutable vector<real> _x;
  /// \brief the sequence of f_x's
  mutable vector<real> _F_x;
  /// \brief the sequence of abs(f_x)'s
  mutable vector<real> _abs_F_x;
  /// \brief the sequence of df_x's
  mutable vector<real> _dF_x;
  /// \brief a boolean which is true, if convergance occurs
  mutable bool _converged;
  /// \brief a boolean which is true if the sequence diverged
  mutable bool _diverged;
  /// \brief true if the final x is on the bounday
  mutable bool _boundary_solution;
  /// \brief the number of times a boundary has been hit
  mutable int _number_boundary_hits;
  /// \brief the number of times the lower boundary has been hit
  mutable bool _lower_boundary_hit;
  /// \brief the number of times the uper boundary has been hit
  mutable bool _upper_boundary_hit;
  /// \brief true if the last iteration hit the boundary
  mutable bool _boundary_hit;

  /// \brief this returns the number of times the bounaries have been hit. If boundaries are hit more than twice, then the loop is stoped.
  int check_boundary( real & x, real df_dx) {

    if (_check_boundary) {
      _lower_boundary_hit = _upper_boundary_hit = false;
      if (x < _min_x) {
	x = _min_x + 2 * _epsilon / df_dx;
	_lower_boundary_hit = true;
	_number_boundary_hits++;
      } else if (x > _max_x) {
	x = _max_x - 2 * _epsilon / df_dx;
	_upper_boundary_hit = true;
	_number_boundary_hits++;
      }
      _boundary_hit = _lower_boundary_hit || _upper_boundary_hit;
    }

    return _number_boundary_hits;
  }

 public:

  /// returns _boundary_solution
  bool get_boundary_solution() const {
    return _boundary_solution;
  }

  /// used to set the second order differential
  void set_d2f( real (functor::*d2f)(real)) {
    _d2f = d2f;
  }

  /// used to turn boundary checking on or off
  void set_check_boundary( bool check_boundary) {
    _check_boundary = check_boundary;
  }

  /// used to set the minimum of the solution set
  void set_min_x( real min_x) {
    _min_x = min_x;
  }

  /// used to set the maximum of the solution set
  void set_max_x( real max_x) {
    _max_x = max_x;
  }

  /// used to retrieve the convergance value
  real get_x_converge() const {
    return _x_converge;
  }

  /// used to retrieve the final value of f
  real get_f_converge() const {
    return _f_converge;
  }

  /// this returns a reference to the parent object with the hidden variables
  functor& get_functor() { return _F;}


  /// This constructor sets everything needed by the algorithm.
  /// @param[in]  X0          This is the starting point of the search. It is stored as the data member, _x0.
  /// @param[in]  Epsilon     The is used to define the end point of the iterative procedure, do_iteration. The loop will terminate  when the abs(F(xi) - fsolve)) is less than Epsilon. It is stored as the data member, _epsilon.
  /// @param[in]  Twice_df        If Twice_df then the derivative of F, at point i, is approximated as (F(x[i])+f(x[i-1])*.5  which approximates the derivative of F over the interval x[i] to x[i-1]. It is stored as the data member, _F.
  /// @param[in]  max_iter        The maximum number of iterations before do_iteration, will terminate. It is stored as the data member _max_iter.
  /// @param[in]  F               This reference to a functor object serves as the object for calling the member function pointers f, df and d2f. It is stored as the data member _F.
  /// @param[in]  f               This pointer to a member function is used to calculate the function to solved. It is stored as the data member _f.
  /// @param[in]  df              This pointer to a member function is used to calculate the derivative of the function to solved. It is stored as the data member _df.
  /// @param[in]  d2f             This pointer to a member function is used to calculate the second derivative of the function to solved. It is stored as the data member _d2f.
  /// @param[in]  fsolve          This is the value of f to be solved. It's default value is 0. It is stored as the data member _fsolve.
  /// @param[in]  check_boundary  If true then the algorithm will check constrain the search so that it only stays within min_x and max_x.  If the search goes outside of the range then the solution is set to the boundary value. It is stored as the data member _check_boundary.
  NewtonRaphsonSolve0(
    real X0,
    real Epsilon,
    bool Twice_df,
    unsigned int max_iter,
    functor& F,
    real (functor::*f)(real),
    real (functor::*df)(real),
    real (functor::*d2f)(real) = NULL,
    real fsolve = 0,
    bool check_boundary = false,
    real min_x = 0,
    real max_x = 0)
    :  _x0( X0), _epsilon( Epsilon), _twice_df( Twice_df), _max_iter(max_iter), _F(F), _f(f), _df(df), _d2f(d2f), _fsolve(fsolve), _check_boundary(check_boundary), _min_x(min_x), _max_x(max_x),
    _x_converge(0), _f_converge(0),
    _x(max_iter), _F_x(max_iter), _abs_F_x(max_iter), _dF_x(max_iter), _converged(false), _diverged(false),
    _lower_boundary_hit(false), _upper_boundary_hit(false), _boundary_hit(false)
  {

  }

  /// This function does the calculations.
  double do_iteration( ostream* output) {
    int i;

    _converged = _diverged = _boundary_solution = false;

    _x[0] = _x0;

    for (i=_number_boundary_hits=0; ;) {
      _F_x[i]       = (_F.*_f)( _x[i]) - _fsolve;

      _abs_F_x[i] = (_F_x[i] > 0)? _F_x[i] : -_F_x[i];

      if (_abs_F_x[i] < _epsilon) {
	_x_converge = _x[i];
	_f_converge = _F_x[i];
	_converged = true;
	break;
      }

      if (i>0 && _abs_F_x[i] > _abs_F_x[i-1]) {
	_diverged = true;
	break;
      }

      if (output) {
	*output << i << " " << _x[i] << " " << _F_x[i] << " " << _dF_x[i] << endl;
      }

      _dF_x[i] = (_F.*_df)( _x[i]);

      ++i;
      if (i >= _max_iter) break;

      if (_d2f != 0) {
	real dx = (_F_x[i-1]) / _dF_x[i-1];
	real df2_dx2 = (_F.*_d2f) ( _x[i-1]);
	_dF_x[i-1] -= .5 * dx * df2_dx2;
      } else if (_twice_df) {
	_x[i] = _x[i-1] - (_F_x[i-1]) / _dF_x[i-1];
	if (check_boundary( _x[i], _dF_x[i-1])>2) {
	  goto boundary;
	}

	if (_boundary_hit) {
	  _dF_x[i-1] = (_F.*_df)( _x[i]);
	} else {
	  real df_dx = (_F.*_df)( _x[i]);
	  _dF_x[i-1] = .5 * (_dF_x[i-1] + df_dx);
	}
      }

      _x[i] = _x[i-1] - (_F_x[i-1]) / _dF_x[i-1];

      if (check_boundary( _x[i], _dF_x[i-1])>2) {
	goto boundary;
      }

      continue;

    boundary:
      if (_lower_boundary_hit) {
	_boundary_solution = true;
	_x_converge = _min_x;
	_f_converge = (_F.*_f)( _min_x);
	break;
      }

      if (_upper_boundary_hit) {
	_boundary_solution = true;
	_x_converge = _max_x;
	_f_converge = (_F.*_f)( _max_x);
	break;
      }
    }

    if (_diverged && output) {
      *output << "answer diverged" << endl;
    }

    if (_boundary_solution && output) {
      *output << i << " " << _x_converge << " " << _f_converge << " boundary_solution" << endl;
        cerr << " boundary_solution" << endl;
        return _x[i];
    } else if (_converged && output) {
      *output << i << " " << _x[i] << " " << _F_x[i] << " answer converged" << endl;
      return _x[i];
    } else if (_converged)
    	return _x[i];
  }
};


#endif
