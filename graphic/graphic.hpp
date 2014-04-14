#ifndef GRAPHIC_HPP
#define GRAPHIC_HPP

#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glut.h"

#include "../bitmap/bitmap.h"
#include "../defines.hpp"

struct Tex {
    GLuint textName;
    int sizX, sizY;
    GLubyte *data;

    void print();
};

class Graphic {
private:
    static std::vector<Tex*> textures;

public:
    static int createTexture(int x, int y);
    static int createTexture(std::string name);
    static void setCamera(const Vec&, const Mat&);
    static void open();
    static void close();
    static void print();
    static void display();
    static void drawRectangle(TYP xmin, TYP xmax, TYP ymin, TYP ymax, int texNum);
};

#endif
