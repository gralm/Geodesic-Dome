#include "controller.hpp"

using namespace std;

Vec *Controller::pos = 0;
Mat *Controller::ori = 0;
int Controller::mouseX = 0, Controller::mouseY = 0;
TYP Controller::rotationVelocity = 0.01;
TYP Controller::distance = 3.500;

bool Controller::privatelyCreatedPosition = false;
bool Controller::privatelyCreatedOrientation = false;
bool Controller::mouseRotationControllerMethod = true;
int Controller::keyControl = 0;

void Controller::motionFunc(int x, int y)
{
//    cout << "Mouse moving at\t[" << x << ", " << y << "]" << endl;

    TYP dx = (x-mouseX)*rotationVelocity;
    TYP dy = (y-mouseY)*rotationVelocity;
    if (mouseRotationControllerMethod) {
        rotateCamera(dx, dy);
    }
    mouseX = x;
    mouseY = y;
    Graphic::setCamera(*pos, *ori);
}

void Controller::mouseFunc(int button, int state, int x, int y)
{
    switch (button) {
        case GLUT_LEFT_BUTTON:
            cout << "Left button \t";
            break;

        case GLUT_MIDDLE_BUTTON:
            cout << "Middle button \t";
            break;

        case GLUT_RIGHT_BUTTON:
            cout << "Right button \t";
            break;
    }

    switch (state) {
        case GLUT_UP:
            cout << "Up";
            break;
        case GLUT_DOWN:
            cout << "Down";
            mouseX = x;
            mouseY = y;
            break;
    }

    cout << "\t at [" << x << ", " << y << "]" << endl;
    //cout << *ori << endl;

}

void Controller::keyboard(unsigned char key, int x, int y)
{
    switch (key) {
        case 27:
        case 'q':
            Graphic::close();
            Controller::close();
            exit(0);
            break;
        case GLUT_KEY_LEFT:
        case GLUT_KEY_DOWN:
        case GLUT_KEY_PAGE_DOWN:
            distance *= 1.0 - rotationVelocity;
            rotateCamera(0.0, 0.0);
            cout << "down" << endl;
            break;
        case GLUT_KEY_RIGHT:
        case GLUT_KEY_PAGE_UP:
        case GLUT_KEY_UP:
            distance *= 1.0 + rotationVelocity;
            rotateCamera(0.0, 0.0);
            cout << "up" << endl;
            break;
        break;
            default:
        break;
    }
    cout << "key: " << key << endl;
}

const Vec& Controller::getCameraPosition()
{
    if (!pos)
        setupPosition();
    return *pos;
}

const Mat& Controller::getCameraOrientation()
{
    if (!ori)
        setupOrientation();
    return *ori;
}

void Controller::setupOrientation()
{
    ori =  new Mat( 0.0, 0., -1.0,
                    0.0, 1.0, 0.0,
                    1.0, 0.0, 0.0);
    privatelyCreatedOrientation = true;
}

void Controller::setupPosition()
{
    pos = new Vec(0., 0., distance);
    privatelyCreatedPosition = true;
}

void Controller::close()
{
    setPosition(0);
    setOrientation(0);
}

void Controller::rotateCamera(TYP dx, TYP dy)
{
    Mat R_(1.0, -dy, dx,
           dy, 1.0, 0.0,
           -dx, 0.0, 1.0);
    *ori = R_ * *ori;
    *pos = -(*ori)[0]*distance;
}

void Controller::setPosition(Vec *pos_)
{
    if (privatelyCreatedPosition){
        delete pos;
        privatelyCreatedPosition = false;
    }
    pos = pos_;
}

void Controller::setOrientation(Mat *ori_)
{
    if (privatelyCreatedOrientation){
        delete ori;
        privatelyCreatedOrientation = false;
    }
    ori = ori_;
}

bool Controller::update()
{
    if (keyControl&CTRL_DOWN) {
        distance *= 1.0 - rotationVelocity;
        *pos = -(*ori)[0]*distance;
        Graphic::setCamera(*pos, *ori);
    }

    if (keyControl&CTRL_UP) {
        distance *= 1.0 + rotationVelocity;
        *pos = -(*ori)[0]*distance;
        Graphic::setCamera(*pos, *ori);
    }
    return false;
}

void Controller::specialFunc(int key, int x, int y)
{
    switch(key) {
        case GLUT_KEY_PAGE_DOWN:
            keyControl ^= CTRL_DOWN;
            break;
        case GLUT_KEY_PAGE_UP:
            keyControl ^= CTRL_UP;
            break;
    }
}

void Controller::specialUpFunc(int key, int x, int y)
{
    switch(key) {
        case GLUT_KEY_PAGE_DOWN:
            keyControl &= ~CTRL_DOWN;
            break;
        case GLUT_KEY_PAGE_UP:
            keyControl &= ~CTRL_UP;
            break;
    }
}


