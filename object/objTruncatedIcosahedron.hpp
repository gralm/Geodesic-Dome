#ifndef __OBJTRUNCATEDICOSAHEDRON_HPP__
#define __OBJTRUNCATEDICOSAHEDRON_HPP__

#include "../defines.hpp"
#include "object.hpp"


class ObjTruncatedIcosahedronFN: public ObjectFN {
public:
	ObjTruncatedIcosahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori);
};

#endif