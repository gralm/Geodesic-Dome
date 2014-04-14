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

#ifndef MATRIX3_H
#define MATRIX3_H

namespace tmath{

template<typename T, int N, int M>
class matrix;

template<typename T>
inline vectorn<T,3>& normalize(vectorn<T,3>& v);

/**
	@author 		Cengiz Terzibas
*/

template<typename T>
class matrix<T,3,3> {
public:
	typedef T value_type;
	T xx, xy, xz;
	T yx, yy, yz;
	T zx, zy, zz;
	// constructors
	matrix<T,3,3>()	:	xx(static_cast<T>(0.0)),
									xy(static_cast<T>(0.0)),
									xz(static_cast<T>(0.0)),
									yx(static_cast<T>(0.0)),
									yy(static_cast<T>(0.0)),
									yz(static_cast<T>(0.0)),
									zx(static_cast<T>(0.0)),
									zy(static_cast<T>(0.0)),
									zz(static_cast<T>(0.0)) {
	}
	matrix<T,3,3>(	T a0, T a1, T a2,
							T a3,	T a4, T a5,
							T a6, T a7, T a8)
						:	xx(a0), xy(a1), xz(a2),
							yx(a3), yy(a4), yz(a5),
							zx(a6), zy(a7), zz(a8) {
	}
	matrix<T,3,3>(const T* m) {
		xx = m[0];	xy = m[1];	xz = m[2];
		yx = m[4];	yy = m[5];	yz = m[6];
		zx = m[8];	zy = m[9];	zz = m[10];
	}
	matrix<T,3,3>(const matrix<T,3,3>& m) {
		xx = m.xx; xy = m.xy; xz = m.xz;
		yx = m.yx; yy = m.yy; yz = m.yz;
		zx = m.zx; zy = m.zy; zz = m.zz;
	}
	// assignment operations
	inline matrix<T,3,3> operator+=(const matrix<T,3,3>& m) {
		xx += m.xx; xy += m.xy; xz += m.xz;
		yx += m.yx; yy += m.yy; yz += m.yz;
		zx += m.zx; zy += m.zy; zz += m.zz;
		return *this;
	}
	inline matrix<T,3,3> operator-=(const matrix<T,3,3>& m) {
		xx -= m.xx; xy -= m.xy; xz -= m.xz;
		yx -= m.yx; yy -= m.yy; yz -= m.yz;
		zx -= m.zx; zy -= m.zy; zz -= m.zz;
		return *this;
	}
	inline matrix<T,3,3> operator*=(const T &num) {
		xx *= num; xy *= num; xz *= num;
		yx *= num; yy *= num; yz *= num;
		zx *= num; zy *= num; zz *= num;
		return *this;
	}
	inline matrix<T,3,3> operator/=(const T& num) {
		xx /= num; xy /= num; xz /= num;
		yx /= num; yy /= num; yz /= num;
		zx /= num; zy /= num; zz /= num;
		return *this;
	}
	// stream operations
	friend std::ostream& operator<<(std::ostream& out, const matrix<T,3,3>& m) {
		out << m.xx << " " << m.xy << " " << m.xz << std::endl;
		out << m.yx << " " << m.yy << " " << m.yz << std::endl;
		out << m.zx << " " << m.zy << " " << m.zz << std::endl;
		return out;
	}
	// comparison operations
	inline const bool operator == (const matrix<T,3,3>& m) const {
		return (	        xx == m.xx && xy == m.xy && xz == m.xz &&
							yx == m.yx && yy == m.yy && yz == m.yz &&
							zx == m.zx && zy == m.zy && zz == m.zz);
	}
	inline const bool operator != ( const matrix<T,3,3>& m) const {
		return !(m == *this);
	}
	// unary operations
	inline const matrix<T,3,3> operator - () const {
		return matrix<T,3,3>(-xx, -xy, -xz,
											-yx, -yy, -yz,
											-zx, -zy, -zz);
	}
	// binary operations
	inline friend matrix<T,3,3> operator+(const matrix<T,3,3>& m1,const matrix<T,3,3>& m2) {
		return matrix<T,3,3>(m1.xx + m2.xx, m1.xy + m2.xy, m1.xz + m2.xz,
											m1.yx + m2.yx, m1.yy + m2.yy, m1.yz + m2.yz,
											m1.zx + m2.zx, m1.zy + m2.zy, m1.zz + m2.zz);
	}
	inline friend matrix<T,3,3> operator-(const matrix<T,3,3>& m1,const matrix<T,3,3>& m2) {
		return matrix<T,3,3>(m1.xx - m2.xx, m1.xy - m2.xy, m1.xz - m2.xz,
											m1.yx - m2.yx, m1.yy - m2.yy, m1.yz - m2.yz,
											m1.zx - m2.zx, m1.zy - m2.zy, m1.zz - m2.zz);
	}
	inline const matrix<T,3,3> operator*(const T& num) const {
		return matrix<T,3,3>(xx * num, xy * num, xz * num,
											yx * num, yy * num, yz * num,
											zx * num, zy * num, zz * num);
	}
	friend inline const matrix<T,3,3> operator*( const T &s, const matrix<T,3,3>& m ) {
		return m * s;
	}
	inline friend matrix<T,3,3> operator*(const matrix<T,3,3>& m1, const matrix<T,3,3>& m2 )  {
		return matrix<T,3,3>(m1.xx * m2.xx + m1.xy * m2.yx + m1.xz * m2.zx,
											m1.xx * m2.xy + m1.xy * m2.yy + m1.xz * m2.zy,
											m1.xx * m2.xz + m1.xy * m2.yz + m1.xz * m2.zz,
											m1.yx * m2.xx + m1.yy * m2.yx + m1.yz * m2.zx,
											m1.yx * m2.xy + m1.yy * m2.yy + m1.yz * m2.zy,
											m1.yx * m2.xz + m1.yy * m2.yz + m1.yz * m2.zz,
											m1.zx * m2.xx + m1.zy * m2.yx + m1.zz * m2.zx,
											m1.zx * m2.xy + m1.zy * m2.yy + m1.zz * m2.zy,
											m1.zx * m2.xz + m1.zy * m2.yz + m1.zz * m2.zz);
	}
	inline friend const vectorn<T,3> operator*(const matrix<T,3,3>& m, const vectorn<T,3>& v) {
		return vectorn<T,3>(v.x*m.xx + v.y*m.xy + v.z*m.xz,
											v.x*m.yx + v.y*m.yy + v.z*m.yz,
											v.x*m.zx + v.y*m.zy + v.z*m.zz);

	}
	inline const matrix<T,3,3> operator/(const T& num) const {
		float val = (static_cast<T>(1.0)/num);
		return matrix<T,3,3>(xx * val, xy * val, xz * val,
											yx * val, yy * val, yz * val,
											zx * val, zy * val, zz * val);
	}
	// indexing operations
	void row(unsigned int pRow, const vectorn<T,3>& pV) {
		unsigned idx = 4*pRow;
		(*this)[idx]		= pV.x;
		(*this)[idx + 1]	= pV.y;
		(*this)[idx + 2]	= pV.z;
	}
	void column(unsigned int pColumn, const vectorn<T,3>& pV) {
		(*this)[pColumn]			= pV.x;
		(*this)[pColumn + 3]	= pV.y;
		(*this)[pColumn + 6]	= pV.z;
	}
	inline const vectorn<T,3> row(int idx) {
		return vectorn<T,3>(*(&xx + 3*idx),*(&xx + 3*idx +1),*(&xx + 3*idx +2));
	}
	inline const vectorn<T,3> column(int idx) {
		return vectorn<T,3>(*(&xx + idx),*(&xx + idx+ 3),*(&xx + idx + 6));
	}
	inline void normalize_column() {
		vectorn<T,3> v0,v1,v2, u0, u1, u2;
		v0 = vectorn<T,3>(xx, yx,zx);
		v1 = vectorn<T,3>(xy, yy,zy);
		v2 = vectorn<T,3>(xz, yz,zz);
		u0 = v0;
		normalize(u0);
		u1 = v1 - ((v1*u0)/(u0*u0))*u0;
		normalize(u1);
		u2 = v2 - ((v2*u0)/(u0*u0))*u0 - ((v2*u1)/(u1*u1))*u1;
		normalize(u2);
		column(0, u0);
		column(1, u0);
		column(2, u0);
	}
	// cast operations
	operator T*() {
		return &xx;
	}
	operator const T*() const {
		return &xx;
	}



        //eget:
	inline friend const matrix<T,3,3> operator&(const vectorn<T,3>& v, const matrix<T,3,3>& m) {
		return matrix<T,3,3>(  v.y*m.xz - v.z*m.xy, v.z*m.xx - v.x*m.xz, v.x*m.xy - v.y*m.xx,
                                v.y*m.yz - v.z*m.yy, v.z*m.yx - v.x*m.yz, v.x*m.yy - v.y*m.yx,
                                v.y*m.zz - v.z*m.zy, v.z*m.zx - v.x*m.zz, v.x*m.zy - v.y*m.zx);
	}
    matrix<T,3,3>(const vectorn<T,3>& a, const vectorn<T,3>& b, const vectorn<T,3>& c) {
		xx = a.x;	xy = a.y;	xz = a.z;
		yx = b.x;	yy = b.y;	yz = b.z;
		zx = c.x;	zy = c.y;	zz = c.z;
	}

        // first order taylor expansion of rotation wih small W, and normalisation and (orthogonalization?)
        // *this is a rotational matrix and will stay so for low all Ws
	inline void quickRotateWC(const vectorn<T,3> &W)
	{
            // X += dt (W x X)
	    T tfx = W.y*xz - W.z*xy;
	    T tfy = W.z*xx - W.x*xz;
	    xz += W.x*xy - W.y*xx;
	    xx += tfx;
	    xy += tfy;

	        // Y += dt (W x Y)
	    tfx = W.y*yz - W.z*yy;
	    tfy = W.z*yx - W.x*yz;
	    yz += W.x*yy - W.y*yx;
	    yx += tfx;
	    yy += tfy;

            // X *= (1.5 - .5X^2)
        tfx = xx*xx + xy*xy + xz*xz;
        tfy = 1.5 - .5*tfx;
        xx *= tfy;  xy *= tfy;  xz *= tfy;

            // Y -= X (X*Y)
        tfx = xx*yx + xy*yy + xz*yz;
        yx -= xx*tfx;   yy -= xy*tfx;   yz -= xz*tfx;

            // Y *= (1.5 - .5X^2)
        tfx = yx*yx + yy*yy + yz*yz;
        tfy = 1.5 - tfx*.5;
        yx *= tfy;  yy *= tfy;  yz *= tfy;

            // Z = X x Y
        zx = xy*yz - xz*yy;
        zy = xz*yx - xx*yz;
        zz = xx*yy - xy*yx;
	}

    inline const vectorn<T,3> operator[] (int n) const {
        switch (n){
            case 0:
                return vectorn<T,3>(xx, xy, xz);
                break;
            case 1:
                return vectorn<T,3>(yx, yy, yz);
                break;
            default:
                return vectorn<T,3>(zx, zy, zz);
                break;
        }
	}

	void set(const int n, const vectorn<T,3>& A) {
	    switch (n){
            case 0:
                xx = A.x;   xy = A.y;   xz = A.z;
                break;
            case 1:
                yx = A.x;   yy = A.y;   yz = A.z;
                break;
            default:
                zx = A.x;   zy = A.y;   zz = A.z;
                break;
        }
	}

    inline const matrix<T,3,3> operator ~() const {
		return matrix<T,3,3>(xx, yx, zx,
                             xy, yy, zy,
                             xz, yz, zz);
    }

    inline void WcRotate(vectorn<T,3>& N_)
    {

        T len = (T) sqrt(N_.x*N_.x + N_.y*N_.y + N_.z*N_.z);
        vectorn<T,3> N(N_/len);

        T s = sin(len);
        T c = cos(len);
        T c1= 1 - c;

        T x, y, f;

            f = c1*(N.x*xx + N.y*xy + N.z*xz);
        x = s*(N.y*xz - N.z*xy) + N.x*f + c*xx;
        y = s*(N.z*xx - N.x*xz) + N.y*f + c*xy;
        xz = s*(N.x*xy - N.y*xx) + N.z*f + c*xz;
        xx = x;     xy = y;

            f = c1*(N.x*yx + N.y*yy + N.z*yz);
        x = s*(N.y*yz - N.z*yy) + N.x*f + c*yx;
        y = s*(N.z*yx - N.x*yz) + N.y*f + c*yy;
        yz = s*(N.x*yy - N.y*yx) + N.z*f + c*yz;
        yx = x;     yy = y;

            f = c1*(N.x*zx + N.y*zy + N.z*zz);
        x = s*(N.y*zz - N.z*zy) + N.x*f + c*zx;
        y = s*(N.z*zx - N.x*zz) + N.y*f + c*zy;
        zz = s*(N.x*zy - N.y*zx) + N.z*f + c*zz;
        zx = x;     zy = y;
    }

    inline void normalize()
    {
        T tf = sqrt(xx*xx + xy*xy + xz*xz);
        if (tf == 0.0) {
            xx = 1.0;
        } else {
            tf = 1.0 / tf;
            xx *= tf;     xy *= tf;   xz *= tf;
        }

        tf = xx*yx + xy*yy + xz*yz;
        if (tf > .999 || tf < -.999){
            if (xx < 0.9 && xx > -0.9){
                yx = 1.0;   yy = 0.0;   yz = 0.0;
            } else {
                yx = 0.0;   yy = 1.0;   yz = 0.0;
            }
            tf = xx*yx + xy*yy + xz*yz;
        }

        yx -= tf*xx;    yy -= tf*xy;    yz -= tf*xz;
        tf = 1.0 / sqrt(yx*yx + yy*yy + yz*yz);
        yx *= tf;   yy *= tf;   yz *= tf;

        zx = xy*yz - xz*yy;
        zy = xz*yx - xx*yz;
        zz = xx*yy - xy*yx;
    }

	inline friend vectorn<T,3> operator* (const vectorn<T,3>& V, const matrix<T,3,3>& M)
	{
	    vectorn<T,3> Vret(  V.x*M.xx + V.y*M.yx + V.z*M.zx,
                            V.x*M.xy + V.y*M.yy + V.z*M.zy,
                            V.x*M.xz + V.y*M.yz + V.z*M.zz);
	    return Vret;
	}

};


}

#endif

