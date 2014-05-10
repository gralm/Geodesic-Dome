#ifndef __OBJCUBE_HPP__
#define __OBJCUBE_HPP__

#include "../defines.hpp"
#include "object.hpp"


class ObjCubeFN: public ObjectFN {
public:
	ObjCubeFN(const Vec &Pos, const Vec &Siz, const Mat &Ori);
};

#endif