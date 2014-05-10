#ifndef __OBJTETRAHEDRON_HPP__
#define __OBJTETRAHEDRON_HPP__

#include "../defines.hpp"
#include "object.hpp"


class ObjTetrahedronFN: public ObjectFN {
public:
	ObjTetrahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori);
};

#endif

