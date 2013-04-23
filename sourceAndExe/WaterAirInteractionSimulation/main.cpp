
#include "FluidSimulation.h"


GLint window_width = 512, window_height = 512;

float *png_data;
unsigned int png_width, png_height;



GLfloat startx, starty;
GLfloat model_angle1 = 0.0, model_angle2 = 0.0, scale = 1.0, eye[3] = { 0.0, 0.5, 1.373+0.5 + 0.1 };
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
bool is_save = false;


void keyboard( unsigned char c, int x, int y )
{
    switch(c){
    case 'q':
        exit(0);
    case 'z':
        K++;
        break;
    case 'x': 
        K--;
        break;
    case 'a':
        Mu_w++;
        break;
    case 's': 
        Mu_w--;
        break;
    case 'r': 
        is_save = true;
        break;
    case 'g': 
        is_proceed = true;
        break;
    }

    glutPostRedisplay();
}



void renderBitmapString(float x, float y, float z, void *font, char *string) 
{  
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glRasterPos3f(x, y, z);

    int len = (int) strlen(string);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(font, string[i]);
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
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
    static float *buffer;

    static bool init = true;
    if(init){
        sph.InitializeGrids();

        buffer = new float [window_width*window_height*3];

        init = false;
    }

    glClearColor( 1.0, 1.0, 1.0, 0.0 );
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2],
              eye[0], eye[1], 0.0,
                 0.0,    1.0,   0.0 );
    glScalef(scale, scale, scale);
    glRotatef(model_angle1, 0, 1, 0);
    glRotatef(model_angle2, 1, 0, 0);
    
    glEnable(GL_LIGHTING);

//    if(is_proceed){
        for(int i = 0; i < 33; i++) sph.ProcessStep();

        sph.RenderParticles();
        sph.OutputParticleLocation();
//    }

    glDisable(GL_LIGHTING);

    double bbmin[5][3] = { {-0.63,  0.0, 0.0}, {-0.45, 0.18, -0.0, },  {-0.09, -0.0, -0.0, }, { 0.27, 0.18, -0.0, }, {-0.27, 0.72, -0.0, } };
    double bbmax[5][3] = { { 0.63, 0.90, 0.18}, {-0.27, 0.90,  0.18, }, { 0.09, 0.36,  0.18, }, { 0.45, 0.90,  0.18, }, { 0.27, 0.90,  0.18, },} ;

    glLineWidth(3);
    glColor3d(0, 0, 1);
    DrawBoundingBox(bbmin[0], bbmax[0]);

    glColor3d(0, 1, 0);
    for(int i = 1; i <= 4; i++){
        if(!is_air_trap && i == 4) break;
        DrawBoundingBox(bbmin[i], bbmax[i]);
    }

    
    char info1[256];
    sprintf_s(info1, sizeof(info1), "K: %lf Mu_w: %lf Mu_a: %lf", K, Mu_w, Mu_a);

    glColor3d(0, 0, 0);
    renderBitmapString(-0.95, 0.9, 0, GLUT_BITMAP_HELVETICA_18,  info1);

    glutSwapBuffers();
    glutPostRedisplay();
 }


int main(int argc, char *argv[])
{
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

