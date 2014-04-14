#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <list>
#include <vector>
#include <iostream>
#include <cmath>
#include "tm/tm.h"

#define UPD_GRAPH       40
#define UPD_PHYS        10
#define UPD_PRINT       1000

    // GENERAL DEFINES
#define TYP     double

#define Vec     tmath::vectorn<TYP, 3>
#define Mat     tmath::matrix<TYP, 3, 3>
#define ONE     tmath::matrix<TYP, 3, 3>(1,0,0, 0,1,0, 0,0,1)


#endif
