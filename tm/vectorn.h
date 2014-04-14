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

#ifndef VECTORN_H
#define VECTORN_H

#include <vector>

namespace tmath{

/**
	@class 			vectorn
	@brief 			Class that represents a n-dimensional vector
	@author 		Cengiz Terzibas
*/

template<typename T, int NUM>
class vectorn {
public:
	typedef T value_type;
	T cmp[NUM];
	// constructors
	vectorn<T,NUM>() {
		for(int idx = 0; idx < NUM; ++idx)
			cmp[idx] = static_cast<T>(0.0);
	}
	vectorn<T,NUM>(const T* vector){
		for(int idx = 0; idx < NUM; ++idx)
			cmp[idx] = vector[idx];
	}
	vectorn<T,NUM>(const vectorn<T,NUM>& v)  {
		for(int idx = 0; idx < NUM; ++idx)
			cmp[idx] = v.cmp[idx];
	}
	// assignment operations
	inline const vectorn<T,NUM> operator = ( const vectorn<T,NUM>& v ) {
		for(int idx = 0; idx < NUM; ++idx)
			cmp[idx] = v.cmp[idx];
		return *this;
	}
	inline const vectorn<T,NUM> operator+=(const vectorn<T,NUM>& v) {
		for(int idx = 0; idx < NUM; ++idx)
			cmp[idx] += v.cmp[idx];
		return *this;
	}
	inline const vectorn<T,NUM> operator-=(const vectorn<T,NUM>& v) {
		for(int idx = 0; idx < NUM; ++idx)
			cmp[idx] -= v.cmp[idx];
		return *this;
	}
	inline const vectorn<T,NUM> operator*=(const T& num) {
		for(int idx = 0; idx < NUM; ++idx)
			cmp[idx] *= num;
		return *this;
	}
	inline const vectorn<T,NUM> operator/=(const T& num) {
		const T r = static_cast<T>(1.0)/num;
		for(int idx = 0; idx < NUM; ++idx)
			cmp[idx] *= r;
		return *this;
	}
	// stream operations
	friend std::ostream& operator<<(std::ostream& out, const vectorn<T,NUM>& v) {
		for(int idx = 0; idx < NUM; ++idx)
			out << v.cmp[idx] << " ";
		return out;
	}
	// comparison operations
	inline const bool operator == ( const vectorn<T,NUM>& v ) const {
		bool res = true;
		for(int idx = 0; idx < NUM; ++idx)
			res &= v.cmp[idx] == cmp[idx];
		return res;
	}
	inline const bool operator != ( const vectorn<T,NUM>& v ) const {
		return !(v == *this);
	}
	// unary operations
	inline const vectorn<T,NUM> operator - () const {
		vectorn<T,NUM> v;
		for(int idx = 0; idx < NUM; ++idx)
			v.cmp[idx] = -cmp[idx];
		return v;
	}
	// binary operations
	inline friend const vectorn<T,NUM> operator+(const vectorn<T,NUM>& v1,const vectorn<T,NUM>& v2) {
		vectorn<T,NUM> res;
		for(int idx = 0; idx < NUM; ++idx)
			res[idx] = v1.cmp[idx] + v2.cmp[idx];
		return res;
	}
	inline friend const vectorn<T,NUM> operator-(const vectorn<T,NUM>& v1, const vectorn<T,NUM>& v2) {
		vectorn<T,NUM> res;
		for(int idx = 0; idx < NUM; ++idx)
			res[idx] = v1.cmp[idx] - v2.cmp[idx];
		return res;
	}
	inline friend const T	operator*(const vectorn<T,NUM>& v1,const vectorn<T,NUM>& v2) {
		T res = static_cast<T>(0.0);
		for(int idx = 0; idx < NUM; ++idx)
			res += v1.cmp[idx] * v2.cmp[idx];
		return res;
	}
	inline const vectorn<T,NUM> operator*(const T& num) const {
		vectorn<T,NUM> v;
		for(int idx = 0; idx < NUM; ++idx)
			v.cmp[idx] = cmp[idx]*num;
		return v;
	}
	friend inline const vectorn<T,NUM> operator * ( const T& s, const vectorn<T,NUM>& v ) {
		return v * s;
	}
	inline const vectorn<T,NUM> operator/(const T& num) const {
		vectorn<T,NUM> v;
		const T r = static_cast<T>(1.0)/num;
		for(int idx = 0; idx < NUM; ++idx)
			v.cmp[idx] = cmp[idx]*r;
		return v;
	}
	size_t size() const {
		return NUM;
	}
	// cast operations
	operator T*() {
		return cmp;
	}
	operator const T*() const {
		return cmp;
	}
};

}

#endif

