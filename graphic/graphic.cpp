#include "graphic.hpp"

using namespace std;

vector<Tex*> Graphic::textures;


void Tex::print()
{
    cout << "Tex-name: " << textName << "\t(" << sizX << "|" << sizY << ")" << endl;
}

void Graphic::setCamera(const Vec &P, const Mat &M)
{
    glClear (GL_COLOR_BUFFER_BIT);
    glColor3f (1.0, 1.0, 1.0);
    glLoadIdentity (); /* clear the matrix */
    /* viewing transformation */
    //      detta kommando funkar inte, glömt nån glufil?
    gluLookAt(P.x, P.y, P.z, P.x+M.xx, P.y+M.xy, P.z+M.xz, M.yx, M.yy, M.yz);
    glScalef (1.0, 1.0, 1.0); /* modeling transformation */
    glutWireCube (1.0);
    glFlush ();
}

void Graphic::display()
{
    static TYP asdf = 0.0;
    asdf += 0.001;

    ObjectFN *ofn = World::getAnObject();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    drawRectangle(asdf, asdf+1.0, 0.0, 1.0, 0);
    drawRectangle(-1.0, 0.0, asdf+0.0, asdf+1.0, 1);
    drawObject(ofn);

    glFlush();
    glutSwapBuffers();
}

void Graphic::drawRectangle(TYP xmin, TYP xmax, TYP ymin, TYP ymax, int texNum)
{

    glBindTexture(GL_TEXTURE_2D, textures[texNum]->textName);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex3f(xmin, ymin, 0.0);
    glTexCoord2f(0.0, 1.0); glVertex3f(xmin, ymax, 0.0);
    glTexCoord2f(1.0, 1.0); glVertex3f(xmax, ymax, 0.0);
    glTexCoord2f(1.0, 0.0); glVertex3f(xmax, ymin, 0.0);
    glEnd();
}

void Graphic::drawObject(ObjectFN *ofn)
{
    //cout << "kom hit" << endl;
    glEnableClientState(GL_NORMAL_ARRAY);
    int numOfFaces = 0;
    Face *Fac_ = ofn->getFaces(numOfFaces);

    //cout << "antal faces: " << numOfFaces << endl;

    for (int f=0; f<numOfFaces; f++)
    {
        Edge *from = Fac_[f].from;
        //glBegin(GL_QUADS); 
        glBegin(GL_POLYGON);
        //glBegin(GL_TRIANGLES);
        glNormal3f(Fac_[f].Norm.x, Fac_[f].Norm.y, Fac_[f].Norm.z);
        do {
            Vertex *VVVVVVV = from->fr;
            Vec *asdf = &VVVVVVV->X;
            Vec *v_ = &from->fr->X;
            glVertex3f(v_->x, v_->y, v_->z);
            from = from->next;
        } while (from != Fac_[f].from);
        glEnd();
    }
    glDisableClientState(GL_NORMAL_ARRAY);

}

void Graphic::close()
{
    vector<Tex*>::iterator itT = textures.begin();
    while (!textures.empty())
    {
        Tex *texToDelete = static_cast<Tex*>(*itT);
        delete texToDelete->data;
        delete texToDelete;
        textures.erase(itT);
    }
}

int Graphic::createTexture(int x, int y)
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);

    Tex *newText = new Tex;
    newText->data = new GLubyte[x*y*4];
    newText->sizX = x;
    newText->sizY = y;

    int i, j, c = 0;
    for (i = 0; i < y; i++) {
        for (j = 0; j < x; j++) {
            newText->data[c*4 + 0] = (GLubyte) ((i+j)*32)%255;
            newText->data[c*4 + 1] = (GLubyte) (j*4)%255;
            newText->data[c*4 + 2] = (GLubyte) (i*4)%255;
            newText->data[c*4 + 3] = (GLubyte) 255;
            c++;
        }
    }


    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &newText->textName);
    glBindTexture(GL_TEXTURE_2D, newText->textName);

        // GL_REPEAT om det ska repeteras i de olika rikningarna, annars GL_CLAMP
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

        // vid  Magnification, minification beror på avståndet till det som ska rendreras
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, x, y, 0, GL_RGBA, GL_UNSIGNED_BYTE, newText->data);

    textures.push_back(newText);
    return textures.size()-1;
}

int Graphic::createTexture(std::string fileName)
{
    Bitmap bmp;
    if (!bmp.loadBMP(fileName.c_str()))
    {
        cout << "no such file: " << fileName << endl;
        return -1;
    }

    Tex *newTex = new Tex;
    newTex->sizX = bmp.width;
    newTex->sizY = bmp.height;

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);

    int siz = newTex->sizX*newTex->sizY;
    newTex->data = new GLubyte[siz*4];

    for (int c=0; c<siz; c++) {
        newTex->data[c*4 + 0] = (GLubyte) bmp.data[c*3 + 0];
        newTex->data[c*4 + 1] = (GLubyte) bmp.data[c*3 + 1];
        newTex->data[c*4 + 2] = (GLubyte) bmp.data[c*3 + 2];
        newTex->data[c*4 + 3] = (GLubyte) 255;
    }

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &newTex->textName);
    glBindTexture(GL_TEXTURE_2D, newTex->textName);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, newTex->sizX, newTex->sizY, 0, GL_RGBA, GL_UNSIGNED_BYTE, newTex->data);

    textures.push_back(newTex);
    return textures.size()-1;
}

void Graphic::print()
{
    for (unsigned int i=0; i<textures.size(); i++)
        textures[i]->print();
}
