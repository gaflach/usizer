#ifndef EDGE_ARRAY_H
#define	EDGE_ARRAY_H

#include <algorithm>
#include "fmath.hpp"
using std::swap;

#define MAKE_SELF_OPERATOR( OP ) \
friend void operator OP ( EdgeArray<T> &v0, const EdgeArray<T> v1 ) { v0[RISE] OP v1[RISE], v0[FALL] OP v1[FALL]; } \
friend void operator OP ( EdgeArray<T> &v0, const T            v1 ) { v0[RISE] OP v1; v0[FALL] OP v1; }

#define MAKE_OPERATOR( OP ) \
friend EdgeArray<T> operator OP ( const EdgeArray<T> v0, const EdgeArray<T> v1 ) { return EdgeArray<T>(v0[RISE] OP v1[RISE], v0[FALL] OP v1[FALL]); } \
friend EdgeArray<T> operator OP ( const T            v0, const EdgeArray<T> v1 ) { return EdgeArray<T>(v0       OP v1[RISE], v0       OP v1[FALL]); } \
friend EdgeArray<T> operator OP ( const EdgeArray<T> v0, const T            v1 ) { return EdgeArray<T>(v0[RISE] OP v1      , v0[FALL] OP v1      ); }

template<typename T>
class EdgeArray {

	friend ostream &operator<<(ostream &out, const EdgeArray<T> array ) {
		return out << "(" << array[RISE] << ", " << array[FALL] << ")";
	} // end operator

	MAKE_OPERATOR(+);
	MAKE_OPERATOR(-);
	MAKE_OPERATOR(*);
	MAKE_OPERATOR(/);

	MAKE_SELF_OPERATOR(+=);
	MAKE_SELF_OPERATOR(-=);
	MAKE_SELF_OPERATOR(*=);
	MAKE_SELF_OPERATOR(/=);

	friend EdgeArray<T> operator-( const EdgeArray<T> &v0 ) { return EdgeArray<T>(-v0[RISE], -v0[FALL]); }

	friend EdgeArray<T> max( const EdgeArray<T> v0, const EdgeArray<T> v1 ) { return EdgeArray<T>(max(v0[RISE],v1[RISE]),max(v0[FALL],v1[FALL])); }
	friend EdgeArray<T> min( const EdgeArray<T> v0, const EdgeArray<T> v1 ) { return EdgeArray<T>(min(v0[RISE],v1[RISE]),min(v0[FALL],v1[FALL])); }

	friend EdgeArray<T> abs( const EdgeArray<T> v ) { return EdgeArray<T>(fabs(v[RISE]),fabs(v[FALL])); }
	friend EdgeArray<T> pow( const EdgeArray<T> v, const double exp ) { return EdgeArray<T>(pow(v[RISE], exp),pow(v[FALL], exp)); }
	friend EdgeArray<T> pow2( const EdgeArray<T> v ) { return EdgeArray<T>(v[RISE]*v[RISE],v[FALL]*v[FALL]); }
	friend EdgeArray<T> pow3( const EdgeArray<T> v ) { return EdgeArray<T>(v[RISE]*v[RISE]*v[RISE],v[FALL]*v[FALL]*v[FALL]); }
	friend EdgeArray<T> sqrt( const EdgeArray<T> v) { return EdgeArray<T>(sqrt(v[RISE]),sqrt(v[FALL])); }
	friend EdgeArray<T> exp( const EdgeArray<T> v) { return EdgeArray<T>(fmath::expd(v[RISE]),fmath::expd(v[FALL])); }

private:
	T clsValue[2];
public:
	T &operator[]( const EdgeType edgeType ) {return clsValue[edgeType];}
	T  operator[]( const EdgeType edgeType ) const {return clsValue[edgeType];}

	EdgeArray &operator=(const EdgeArray &array) {
		clsValue[RISE] = array[RISE];
		clsValue[FALL] = array[FALL];
		return *this;
	} // end operator

	EdgeArray(const T rise, const T fall ) {
		clsValue[RISE] = rise;
		clsValue[FALL] = fall;
	} // end constructor

	EdgeArray(){};

	void set( const T rise, const T fall ) {
		clsValue[RISE] = rise;
		clsValue[FALL] = fall;
	} // end method

	T getMax() const { return max(clsValue[RISE], clsValue[FALL]); }
	T getMin() const { return min(clsValue[RISE], clsValue[FALL]); }
	T getRise() const { return clsValue[RISE]; }
	T getFall() const { return clsValue[FALL]; }

	EdgeArray<T> getReversed() const { return EdgeArray(getFall(), getRise()); }

	void reverse() { return swap(clsValue[RISE], clsValue[FALL]); }

	double aggregate() const { return clsValue[RISE] + clsValue[FALL]; }
}; // end class

#endif

