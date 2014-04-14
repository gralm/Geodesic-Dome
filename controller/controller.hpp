#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "../defines.hpp"
#include "../graphic/graphic.hpp"

#define CTRL_DOWN       1
#define CTRL_UP         2

class Controller {
private:
    static Vec *pos;
    static Mat *ori;

    static int mouseX, mouseY;
    static TYP rotationVelocity;
    static TYP distance;
    static int keyControl;

    static bool privatelyCreatedPosition;
    static bool privatelyCreatedOrientation;
    static bool mouseRotationControllerMethod;

    static void setupOrientation();
    static void setupPosition();
    static void rotateCamera(TYP dx, TYP dy);
    static void close();

public:
    static void mouseFunc(int button, int state, int x, int y);
    static void motionFunc(int x, int y);
    static void keyboard (unsigned char key, int x, int y);
    static void specialFunc(int key, int x, int y);
    static void specialUpFunc(int key, int x, int y);

    static const Vec& getCameraPosition();
    static void setPosition(Vec*);
    static const Mat& getCameraOrientation();
    static void setOrientation(Mat*);

    static bool update();
};

#endif

