/*
	Copyright (c) 2003-2005 Cengiz Terzibas

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Cengiz Terzibas         cengiz@terzibas.de
*/

#ifndef VECTOR2_H
#define VECTOR2_H

namespace tmath{

template<typename T, int NUM>
class vectorn;

/**
	@author 		Cengiz Terzibas
*/

template<typename T>
class vectorn<T,2> {
public:
	typedef T value_type;
	T x, y;
	// constructors
	vectorn<T,2>() :	x(static_cast<T>(0.0)),
									y(static_cast<T>(0.0)) {
	}
	vectorn<T,2>(T x_, T y_) : x(x_), y(y_) {
	}
	vectorn<T,2>(const T* v) : x(v[0]), y(v[1]) {
	}
	vectorn<T,2>(const vectorn<T,2>& v) : x(v.x), y(v.y) {
	}
	// assignment operations
	inline const vectorn<T,2> operator+=(const vectorn<T,2>& v) {
		x += v.x; y += v.y; return *this;
	}
	inline const vectorn<T,2> operator-=(const vectorn<T,2>& vVector) {
		x -= vVector.x; y -= vVector.y; return *this;
	}
	inline const vectorn<T,2> operator*=(const T& num) {
		x *= num; y *= num;  return *this;
	}
	inline const vectorn<T,2> operator/=(const T& num) {
		const T r = 1.0/num; x *= r; y *= r; return *this;
	}
	// stream operations
	inline friend std::ostream& operator<<(std::ostream& out, const vectorn<T,2>& v) {
		out << v.x << " " << v.y;  return out;
	}
	// comparison operations
	inline const bool operator == ( const vectorn<T,2>& v ) const {
		return (v.x==x && v.y==y );
	}
	inline const bool operator != ( const vectorn<T,2>& v ) const {
		return !(v == *this);
	}
	// unary operations
	inline const vectorn<T,2> operator - () const {
		return vectorn<T,2>( -x, -y);
	}
	// binary operations
	inline friend const vectorn<T,2> operator+(const vectorn<T,2>& v1, const vectorn<T,2>& v2){
		return vectorn<T,2>(v1.x + v2.x, v1.y + v2.y);
	}
	inline friend const vectorn<T,2> operator-(const vectorn<T,2>& v1, const vectorn<T,2>& v2) {
		return vectorn<T,2>(v1.x - v2.x, v1.y - v2.y);
	}
	inline friend const T	operator*(const vectorn<T,2>& v1,const vectorn<T,2>& v2) {
		return T(v1.x * v2.x + v1.y *v2.y );
	}
	inline const vectorn<T,2> operator*(const T& num) const {
		return vectorn<T,2>(x * num, y * num);
	}
	friend inline const vectorn<T,2> operator * ( const T& s, const vectorn<T,2>& v ) {
		return v * s;
	}
	inline const vectorn<T,2> operator/(const T& num) const {
		const T r = 1.0/num; return vectorn<T,2>(x * r, y * r);
	}
	size_t size() const {
		return 2;
	}
	// cast operations
	operator T*() {
		return &x;
	}
	operator const T*() const	{
		return &x;
	}
};

}

#endif

