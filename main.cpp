#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>

#include "graphic/graphic.hpp"
#include "controller/controller.hpp"

class Graphic;

void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum (-1.0, 1.0, -1.0, 1.0, 1.5, 20.0);
    //gluPerspective(60.0, (GLfloat) w/(GLfloat) h, 1.0, 30.0);
    glMatrixMode(GL_MODELVIEW);
    /*glLoadIdentity();
    glTranslatef(0.0, 0.0, -3.6);*/
}

void idle()
{
    static int var = -1;                     //Warning
        // Den här skiten måste köras här MITT I main-loopen för att programmet inte har någon
        //  uppfattning om vilken upplösning som används före glutMainLoop har startats. FAN!!!
    if (var==-1)
    {
        var = 0;
        /*TYP Pf[4];
        int Pi[4];
        glGetIntegerv(GL_VIEWPORT, Pi);
        Graphic::PosPixToFloat(0, 0, Pf);
        Graphic::PosPixToFloat((int) Pi[2], (int) Pi[3], Pf+2);*/
    }


    int time = glutGet(GLUT_ELAPSED_TIME);
    static int time_graph = time;
    static int time_phys = time;
    static int time_print = time;

    if (time > time_graph + UPD_GRAPH)
    {
        Controller::update();
        glutPostRedisplay();
        time_graph += UPD_GRAPH;
        if (var&1){
            printf("too high graph-update\n");
            time_graph = time;
        }
        var |= 1;
    }else
        var &= ~1;


    if (time > time_phys + UPD_PHYS)
    {
        time_phys += UPD_PHYS;
        //Graphic::Update(UPD_PHYS * .001);

        if (var&2){
            printf("too high phys-update\n");
            time_phys = time;
        }
        var |= 2;
    }else
        var &= ~2;


    if ((time > time_print + UPD_PRINT) && UPD_PRINT)
    {
        time_print += UPD_PRINT;
        //Group::Skriv();

        if (var&4){
            printf("too high print-update\n");
            time_print = time;
        }
        var |= 4;
    }else
        var &= ~4;
}


int main(int argc, char** argv)
{
    /*Vec pos = Vec(0., 0., 2.5);
    Mat ori = Mat(0., 0., -1,
                  0, 1.0, 0.,
                  1.0, .0, .0);*/

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(250, 250);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(argv[0]);

    Graphic::setCamera(Controller::getCameraPosition(), Controller::getCameraOrientation());
    Graphic::createTexture(64, 64);
    Graphic::createTexture("pics/meny.bmp");

    glEnable(GL_TEXTURE_2D);
    glutDisplayFunc(Graphic::display);
    glutIdleFunc(idle);
    //glutSetKeyRepeat(GLUT_KEY_REPEAT_OFF);
    glutIgnoreKeyRepeat(true);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(Controller::keyboard);
    glutMouseFunc(Controller::mouseFunc);
    glutMotionFunc(Controller::motionFunc);
    glutSpecialFunc(Controller::specialFunc);
    glutSpecialUpFunc(Controller::specialUpFunc);

    glutMainLoop();
    return 0;
}
