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

#ifndef VECTOR4_H
#define VECTOR4_H

namespace tmath{

template<typename T>
class quaternion;

template<typename T, int NUM>
class vectorn;

/**
	@author 		Cengiz Terzibas
*/

template<typename T>
class vectorn<T,4> {
public:
	typedef T value_type;
	T x, y, z, w;
	// constructors
	vectorn<T,4>() :	x(static_cast<T>(0.0)),
									y(static_cast<T>(0.0)),
									z(static_cast<T>(0.0)),
									w(static_cast<T>(0.0)) {
	}
	vectorn<T,4>(T vx, T vy, T vz, T vw) : x(vx), y(vy), z(vz), w(vw) {
	}
	vectorn<T,4>(const T* v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {
	}
	vectorn<T,4>(const vectorn<T,4>& v) : x(v.x), y(v.y), z(v.z), w(v.w) {
	}
	// assignment operations
	inline const vectorn<T,4> operator+=(const vectorn<T,4>& v) {
		x += v.x; y += v.y; z += v.z; w += v.w; return *this;
	}
	inline const vectorn<T,4> operator-=(const vectorn<T,4>& v) {
		x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this;
	}
	inline const vectorn<T,4> operator*=(const T& num) {
		x *= num; y *= num; z *= num; w *= num; return *this;
	}
	inline const vectorn<T,4> operator/=(const T& num) {
		const T r = (T)1.0/num; x *= r; y *= r; z *=r; w *= r; return *this;
	}
	// stream operations
	inline friend std::ostream& operator<<(std::ostream& out, const vectorn<T,4>& v) {
		out << v.x << " " << v.y << " " << v.z << " " << v.w; return out;
	}
	// comparison operations
	inline const bool operator == ( const vectorn<T,4>& v ) const {
		return (v.x==x && v.y==y && v.z==z && v.w==w);
	}
	inline const bool operator != ( const vectorn<T,4>& v ) const {
		return !(v == *this);
	}
	// unary operations
	inline const vectorn<T,4> operator - () const {
		return vectorn<T,4>( -x, -y, -z, -w);
	}
	// binary operations
	inline friend const vectorn<T,4> operator+(const vectorn<T,4>& v1, const vectorn<T,4>& v2)  {
		return vectorn<T,4>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z,v1.w + v2.w);
	}
	inline friend const vectorn<T,4> operator-(const vectorn<T,4>& v1, const vectorn<T,4>& v2) {
		return vectorn<T,4>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w);
	}
	inline friend const T	operator*(const vectorn<T,4>& v1,const vectorn<T,4>& v2) {
		return T(v1.x * v2.x + v1.y *v2.y + v1.z * v2.z + v1.w * v2.w);
	}
	inline const vectorn<T,4> operator*(const T& num) const {
		return vectorn<T,4>(x * num, y * num, z * num, w * num);
	}
	friend inline const vectorn<T,4> operator * ( const T& s, const vectorn<T,4>& v ) {
		return v * s;
	}
	inline const vectorn<T,4> operator/(const T& num) const {
		const T r = (T)1.0/num; return vectorn<T,4>(x * r, y * r, z * r, w * r);
	}
	size_t size() const {
		return 4;
	}
	// cast operations
	operator T*() {
		return &x;
	}
	operator const T*() const {
		return &x;
	}
};

}

#endif

