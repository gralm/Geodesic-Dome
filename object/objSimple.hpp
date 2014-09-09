#ifndef __OBJSIMPLE_HPP__
#define __OBJSIMPLE_HPP__

#include "../defines.hpp"
#include "object.hpp"


class ObjSimpleFN: public ObjectFN {
public:
	ObjSimpleFN(const Vec &Pos, const Vec &Siz, const Mat &Ori);
};

#endif