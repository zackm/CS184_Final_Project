
#include "FluidSimulation.h"


GLint window_width = 512, window_height = 512;

float *png_data;
unsigned int png_width, png_height;


GLfloat startx, starty;
GLfloat model_angle1 = 0.0, model_angle2 = 0.0, scale = 1.0, eye[3] = { 0.0, 0.5, 1.373+0.5 };
bool    left_click = 0, right_click = 0;


void mouse(int button, int state, int x, int y)
{
    if ( button == GLUT_LEFT_BUTTON ) {
        if ( state == GLUT_DOWN ) {
            left_click = true;
            startx   = x;
            starty   = y;
        }else if (state == GLUT_UP) {
            left_click = false;
        }
    }else{ // button == GLUT_RIGHT_BUTTON
        if ( state == GLUT_DOWN ) {  
            right_click = true;
            startx   = x;
            starty   = y;
        }else if (state == GLUT_UP) {
            right_click = false;
        }
    }

}

void motion( int x, int y )
{ 
    if ( left_click && !right_click ) {       // rotating image
        model_angle1 += (x - startx);
        model_angle2 += (y - starty);
    }else if( !left_click && right_click ) { 
        eye[0] -= (x - startx) / (window_width *0.25);
        eye[1] += (y - starty) / (window_height*0.25);
    }else{ // if( left_click && right_click ) // scaling
        eye[2] -= (y - starty) * 0.01;
    //        scale -= (y - starty) * 0.01;
    }

    startx = x;
    starty = y;

    glutPostRedisplay();
}

bool is_proceed = false;
bool is_render = false;


void keyboard( unsigned char c, int x, int y )
{
    switch(c){
    case 'q':
        exit(0);
    case 'r': 
        if(is_render) is_render = false;
        else          is_render = true;
        break;
    case 'g': 
        is_proceed = true;
        break;
    }

    glutPostRedisplay();
}

void DrawBoundingBox(double *bbmin, double *bbmax)
{
    glPushMatrix();
    glTranslatef((bbmin[0]+bbmax[0])*0.5, (bbmin[1]+bbmax[1])*0.5, (bbmin[2]+bbmax[2])*0.5);
    glScalef(bbmax[0]-bbmin[0], bbmax[1]-bbmin[1], bbmax[2]-bbmin[2]);
    glutWireCube(1.0); 
    glPopMatrix();
}

SPH sph;

void display()
{
    glClearColor( 1.0, 1.0, 1.0, 0.0 );
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2],
              eye[0], eye[1],    0.0,
                 0.0,    1.0,    0.0 );
  //  glScalef(scale, scale, scale);
    glRotatef(model_angle1, 0, 1, 0);
    glRotatef(model_angle2, 1, 0, 0);

    double bbmin[3] = {-0.63,  0.0, 0.0}, bbmax[3] = { 0.63, 0.90, 0.18};

    glLineWidth(3);
    glColor3d(0, 0, 1);
    DrawBoundingBox(bbmin, bbmax);
    
    glPushMatrix();
    glLoadIdentity();
    glRotatef(model_angle1, 0, 1, 0);
    glRotatef(model_angle2, 1, 0, 0);

    double org[3] = { eye[0], eye[1], eye[2], };

    if(is_render){
        sph.InitializeGrids();
        sph.RenderUsingRaycasting(org, window_width, window_height);
    }

    glPopMatrix();

    glutSwapBuffers();
    glutPostRedisplay();
 }


int main(int argc, char *argv[])
{
    GdiplusStartupInput gdiplusStartupInput;
    ULONG_PTR           gdiplusToken;

    // Initialize GDI+.
    GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

    png_data = LoadPngData(L"chiyo_koi2.png", png_width, png_height);

 //   GdiplusShutdown(gdiplusToken);


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(300, 20);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow(argv[0]);
    glutDisplayFunc(display);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    GLInit();
    glutMainLoop();

    return 1;
}

