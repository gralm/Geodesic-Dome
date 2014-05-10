#ifndef __OBJDODECAHEDRON_HPP__
#define __OBJDODECAHEDRON_HPP__

#include "../defines.hpp"
#include "object.hpp"


class ObjDodecahedronFN: public ObjectFN {
public:
	ObjDodecahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori);
};

#endif