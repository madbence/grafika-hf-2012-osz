#include <math.h>
#include <stdlib.h>
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Vector {
   float x, y, z;
   Vector( ) {x = y = z = 0;}
   Vector(float x0, float y0, float z0 = 0) {x = x0; y = y0; z = z0;}
   Vector operator*(float a) {return Vector(x * a, y * a, z * a);}
   Vector operator+(const Vector& v) {return Vector(x + v.x, y + v.y, z + v.z);}
   Vector operator-(const Vector& v) {return Vector(x - v.x, y - v.y, z - v.z);}
   float operator*(const Vector& v) {return (x * v.x + y * v.y + z * v.z);}
   Vector operator%(const Vector& v) {return Vector(y*v.z-z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);}
   float Length() { return sqrt(x * x + y * y + z * z); }
   Vector unit(){return (*this)*(1/Length());}
};
struct Color {
   float r, g, b;
   Color( ) {r = g = b = 0;}
   Color(float r0, float g0, float b0) {r = r0; g = g0; b = b0;}
   Color operator*(float a) {return Color(r * a, g * a, b * a);}
   Color operator*(const Color& c) {return Color(r * c.r, g * c.g, b * c.b);}
   Color operator+(const Color& c) {return Color(r + c.r, g + c.g, b + c.b);}
};
struct Q{
    float s;
    Vector w;
    Q(){};
    Q(float a, const Vector& w):s(a),w(w){};
    Q(float a, float b, float c, float d):s(a),w(b,c,d){};
    Q operator*(Q& q){return Q(s*q.s-w*q.w,q.w*s+w*q.s+w%q.w);}
};
void glCube(){
    glBegin(GL_TRIANGLES);
    int p[8][3]={{-1,-1,-1},{-1,1,-1},{1,1,-1},{1,-1,-1},{-1,-1,1},{-1,1,1},{1,1,1},{1,-1,1}};
    int n[6][3]={{0,0,-1},{0,-1,0},{-1,0,0},{0,1,0},{1,0,0},{0,0,1}};
    int c[]={0,1,2,0,2,3,0,4,7,0,7,3,1,5,4,1,4,0,2,6,5,2,5,1,3,7,6,3,6,2,4,5,6,4,6,7};
    for(int i=0;i<36;i++){
        glNormal3f(n[i/6][0],n[i/6][1],n[i/6][2]);
        glVertex3f(p[c[i]][0],p[c[i]][1],p[c[i]][2]);
    }
    glEnd();
}
Vector sp(float u, float v){
    return Vector(sin(u)*cos(v), cos(u)*cos(v), sin(v));
}
void glSphere(int us=30, int vs=30){
    float du=2*3.1415/us;
    float dv=3.1415/vs;
    glBegin(GL_TRIANGLES);
    for(int uc=0;uc<us;uc++){
        float u0=uc*du,u1=(uc+1)*du;
        for(int vc=0;vc<vs;vc++){
            float v0=-3.1415/2+vc*dv,v1=-3.1415/2+(vc+1)*dv;
            Vector p0=sp(u0, v0),p1=sp(u1, v0),p2=sp(u0, v1),p3=sp(u1, v1);
            glNormal3f(p0.x, p0.y, p0.z);glVertex3f(p0.x, p0.y, p0.z);
            glNormal3f(p2.x, p2.y, p2.z);glVertex3f(p2.x, p2.y, p2.z);
            glNormal3f(p3.x, p3.y, p3.z);glVertex3f(p3.x, p3.y, p3.z);
            glNormal3f(p0.x, p0.y, p0.z);glVertex3f(p0.x, p0.y, p0.z);
            glNormal3f(p3.x, p3.y, p3.z);glVertex3f(p3.x, p3.y, p3.z);
            glNormal3f(p1.x, p1.y, p1.z);glVertex3f(p1.x, p1.y, p1.z);
        }
    }
    glEnd();
}
Vector cp(float u, float v){return Vector(sin(u), cos(u), v*2-1);}
void glCylinder(int su=20, int sv=10){
    float du=3.1415*2/su,dv=1./sv;
    glBegin(GL_TRIANGLES);
    for(int i=0;i<su;i++){
        float u0=i*du,u1=(i+1)*du;
        for(int j=0;j<sv;j++){
            float v0=j*dv,v1=(j+1)*dv;
            Vector p0=cp(u0,v0),p1=cp(u1,v0),p2=cp(u0,v1),p3=cp(u1,v1);
            glNormal3f(p0.x,p0.y,0);glVertex3f(p0.x,p0.y,p0.z);
            glNormal3f(p1.x,p1.y,0);glVertex3f(p1.x,p1.y,p1.z);
            glNormal3f(p3.x,p3.y,0);glVertex3f(p3.x,p3.y,p3.z);
            glNormal3f(p0.x,p0.y,0);glVertex3f(p0.x,p0.y,p0.z);
            glNormal3f(p2.x,p2.y,0);glVertex3f(p2.x,p2.y,p2.z);
            glNormal3f(p3.x,p3.y,0);glVertex3f(p3.x,p3.y,p3.z);
        }
        glNormal3f(0,0,-1);glVertex3f(sin(u0),cos(u0),-1);glVertex3f(sin(u1),cos(u1),-1);glVertex3f(0,0,-1);
        glNormal3f(0,0,1);glVertex3f(sin(u0),cos(u0),1);glVertex3f(sin(u1),cos(u1),1);glVertex3f(0,0,1);
    }
    glEnd();
}
Vector cop(float u, float v){return Vector(sin(u)*(1-v),cos(u)*(1-v),v);}
void glCone(int su=20, int sv=10){
    float du=3.1415*2/su,dv=1./sv;
    glBegin(GL_TRIANGLES);
    for(int i=0;i<su;i++){
        float u0=i*du,u1=(i+1)*du;
        for(int j=0;j<sv;j++){
            float v0=j*dv,v1=(j+1)*dv;
            Vector p0=cop(u0,v0),p1=cop(u1,v0),p2=cop(u0,v1),p3=cop(u1,v1);
            glNormal3f(p0.x,p0.y,(1-v0));glVertex3f(p0.x,p0.y,p0.z);
            glNormal3f(p1.x,p1.y,(1-v0));glVertex3f(p1.x,p1.y,p1.z);
            glNormal3f(p3.x,p3.y,(1-v1));glVertex3f(p3.x,p3.y,p3.z);
            glNormal3f(p0.x,p0.y,(1-v0));glVertex3f(p0.x,p0.y,p0.z);
            glNormal3f(p2.x,p2.y,(1-v1));glVertex3f(p2.x,p2.y,p2.z);
            glNormal3f(p3.x,p3.y,(1-v1));glVertex3f(p3.x,p3.y,p3.z);
        }
        glNormal3f(0,0,-1);glVertex3f(sin(u0),cos(u0),0);glVertex3f(sin(u1),cos(u1),0);glVertex3f(0,0,0);
    }
    glEnd();
}
unsigned int id;
unsigned char tex[256*256*3];
float khaki[]={195./255,176./255,145./255, 1.};
float grey[]={.5,.5,.5,1};
float glass[]={.8,.8,.8,.6};
float spec[]={1,1,1,1};
float red[]={1,0,0,1};
float black[]={0,0,0,1};
const int hs=300;
Vector houses[hs];
void glGenerateAwesomeTexture(unsigned char *t){
    for(int i=0;i<256*256;i++){
        t[i*3]=rand()%128;
        t[i*3+1]=rand()%256;
        t[i*3+2]=rand()%64;
    }
}
void glHouse(){
    glPushMatrix();
    glTranslatef(0,0,1);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,black);
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,grey);
    glCube();
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,red);
    glBegin(GL_TRIANGLES);
        glNormal3f(0,-3,1);glVertex3f(-1,-1,1);glVertex3f(1,-1,1);glVertex3f(0,0,2);
        glNormal3f(3,0,1);glVertex3f(1,-1,1);glVertex3f(1,1,1);glVertex3f(0,0,2);
        glNormal3f(0,3,1);glVertex3f(-1,-1,1);glVertex3f(-1,1,1);glVertex3f(0,0,2);
        glNormal3f(-3,0,1);glVertex3f(-1,1,1);glVertex3f(1,1,1);glVertex3f(0,0,2);
    glEnd();
    glPopMatrix();
}
void glLandscape(int su=50, int sv=50){
    glNormal3f(0,0,1);
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_TRIANGLES);
    for(int i=1;i<=su;i++){
        for(int j=1;j<=sv;j++){
            float u0=1./i,u1=1./(i+1),v0=1./j,v1=1./(j+1);
            glTexCoord2f(u0,v0);glVertex3f(u0-0.5,v0-0.5,0);
            glTexCoord2f(u1,v0);glVertex3f(u1-0.5,v0-0.5,0);
            glTexCoord2f(u0,v1);glVertex3f(u0-0.5,v1-0.5,0);
            glTexCoord2f(u1,v1);glVertex3f(u1-0.5,v1-0.5,0);
            glTexCoord2f(u0,v1);glVertex3f(u0-0.5,v1-0.5,0);
            glTexCoord2f(u1,v0);glVertex3f(u1-0.5,v0-0.5,0);
        }
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
    for(int i=0;i<hs;i++){
        glPushMatrix();
        glTranslatef(houses[i].x,houses[i].y,0);
        glRotatef(houses[i].z*45,0,0,1);
        glScalef(0.005,0.005,(houses[i].z+houses[i].x+houses[i].y)*0.005+0.005);
        glHouse();
        glPopMatrix();
    }
}
void glArrow(){
    glPushMatrix();
        glCylinder();
        glTranslatef(0,0,1);
        glScalef(2,2,1);
        glCone();
    glPopMatrix();
}
void glRotorBlade(){
    glPushMatrix();
        glTranslatef(-2,0,0);
        glScalef(2,0.2,0.05);
        glCube();
    glPopMatrix();
}
void glRotor(int blades=2, float pitch=10){
    glPushMatrix();
        glTranslatef(0,0,-0.5);
        glScalef(0.2,0.2,0.5);
        glCylinder(10,10);
    glPopMatrix();glPushMatrix();
    float da=360./blades;
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, grey);
    for(int i=0;i<blades;i++){
        glPushMatrix();
            glRotatef(pitch, 1,0,0);
            glRotorBlade();
        glPopMatrix();
        glRotatef(da,0,0,1);
    }
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, khaki);
    glPopMatrix();
}
void glLandingSkidPart(){
    glPushMatrix();
        glTranslatef(0,0.7,0);
        glRotatef(90,0,1,0);
        glScalef(0.1,0.1,1.6);
        glCylinder(8,3);
    glPopMatrix();glPushMatrix();
        glTranslatef(0.8,0.7,0);
        glRotatef(20,1,0,0);
        glScalef(0.07,0.07,0.5);
        glTranslatef(0,0,1);
        glCylinder(8,3);
    glPopMatrix();glPushMatrix();
        glTranslatef(-0.8,0.7,0);
        glRotatef(20,1,0,0);
        glScalef(0.07,0.07,0.5);
        glTranslatef(0,0,1);
        glCylinder(8,3);
    glPopMatrix();
}
void glLandingSkids(){
    glPushMatrix();
        glLandingSkidPart();
        glRotatef(180,0,0,1);
        glLandingSkidPart();
    glPopMatrix();
}
void glApache(){
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, khaki);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
    glPushMatrix();
        glScalef(2,0.9,1);
        glSphere(50,50);
    glPopMatrix();glPushMatrix();
        glTranslatef(0,0,-1.2);
        glLandingSkids();
    glPopMatrix();glPushMatrix();
        glRotatef(-10,0,1,0);
        glTranslatef(2,0,0);
        glScalef(2,0.2,0.2);
        glRotatef(90,0,1,0);
        glCylinder();
    glPopMatrix();glPushMatrix();
        glTranslatef(0,0,1.8);
        glRotatef(20,0,0,1);
        glRotor();
    glPopMatrix();glPushMatrix();
        glScalef(0.2,0.2,0.2);
        glRotatef(90,1,0,0);
        glTranslatef(18,3.5,2);
        glRotor(3);
    glPopMatrix();glPushMatrix();
        glTranslatef(3.8,0,1);
        glScalef(0.1,0.1,0.3);
        glRotatef(10,0,1,0);
        glCylinder(10,2);
    glPopMatrix();glPushMatrix();
        glEnable(GL_BLEND);
        glEnable(GL_CULL_FACE);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, glass);
        glTranslatef(-1.2,0,0.6);
        glRotatef(-30,0,1,0);
        glScalef(0.8,0.6,0.4);
        glSphere();
        glDisable(GL_CULL_FACE);
        glDisable(GL_BLEND);
    glPopMatrix();
}
const int screenWidth = 600;
const int screenHeight = 600;
float mrand(){ return (float)rand()/RAND_MAX-0.5;}
void onInitialization( ) {
    glViewport(0, 0, screenWidth, screenHeight);
    glMatrixMode(GL_PROJECTION);
    gluPerspective(60,1,1,1000);
    glMatrixMode(GL_MODELVIEW);
    gluLookAt(-15,0,8,0,0,0,0,0,1);
    glEnable(GL_NORMALIZE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    float ambient[]={.3,.3,.3},diffuse[]={.8,.8,.8,1},specular[]={1,1,1},position[]={-1,1,1,0};
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT0, GL_POSITION, position);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,ambient);
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,100);
    glGenerateAwesomeTexture(tex);
    glGenTextures(1,&id);
    glBindTexture(GL_TEXTURE_2D, id);
    glTexImage2D(GL_TEXTURE_2D,0, GL_RGB, 256, 256, 0,GL_RGB, GL_UNSIGNED_BYTE, tex);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);
    for(int i=0;i<hs;i++){ houses[i]=Vector(mrand(),mrand(),mrand()); }
}
Q q(1,0,0,0);
void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
        glTranslatef(0,0,-10);
        glScalef(200,200,200);
        glLandscape();
    glPopMatrix();glPushMatrix();
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,red);
        Vector w=q.w.unit(),r((w.x+1)/2,w.y/2,w.z/2);
        if(r.Length()<0.01) r=Vector(0,0,1);
        if(q.w.Length()>0.01){
            glRotatef(180,r.x,r.y,r.z);
            glScalef(fabs(q.s)*4,0.15,0.15);
            glTranslatef(1,0,0);
            glRotatef(90,0,1,0);
            glArrow();}
    glPopMatrix();glPushMatrix();
        glRotatef(acos(q.s)/3.1415*360,q.w.x,q.w.y,q.w.z);
        glApache();
    glPopMatrix();
    glutSwapBuffers();
}
void onKeyboard(unsigned char key, int x, int y) {
    float r=10*3.1415/180;
    Q n;
    switch(key){
        case 'r':case 'R':n=Q(cos(r),sin(r),0,0); break;
        case 'e':case 'E':n=Q(cos(-r),sin(-r),0,0); break;
        case 'p':case 'P':n=Q(cos(r),0,sin(r),0); break;
        case 'o':case 'O':n=Q(cos(-r),0,sin(-r),0); break;
        case 'x':case 'X':n=Q(cos(-r),0,0,sin(-r)); break;
        case 'y':case 'Y':n=Q(cos(r),0,0,sin(r)); break;
        default: return;}
    q=q*n;
    glutPostRedisplay();
}
void onMouse(int button, int state, int x, int y) {}
void onIdle( ) {}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Grafika funky feladat");
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    onInitialization();
    glutDisplayFunc(onDisplay);
    glutMouseFunc(onMouse);
    glutIdleFunc(onIdle);
    glutKeyboardFunc(onKeyboard);
    glutMainLoop();
    return 0;
}