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
   Vector( ) { x = y = z = 0; }
   Vector(float x0, float y0, float z0 = 0) { x = x0; y = y0; z = z0; }
   Vector operator*(float a) { return Vector(x * a, y * a, z * a); }
   Vector operator+(const Vector& v) { return Vector(x + v.x, y + v.y, z + v.z); }
   Vector operator-(const Vector& v) { return Vector(x - v.x, y - v.y, z - v.z); }
   float operator*(const Vector& v) { return (x * v.x + y * v.y + z * v.z); }
   Vector operator/(float a) { return Vector(x/a, y/a, z/a); }};
class Spline{
    protected:
    const static int pMax=100, res=1000;
    const Vector lb, rt;
    Vector points[pMax];
    double time[pMax];
    int size;
    public:
    Spline():size(0),lb(0,0),rt(1000,1000){};
    void redistribute(){
        double t=0,dt=0,ddt=2./(size-1.)/size;
        for(int i=0;i<size;i++,dt+=ddt,t+=dt) time[i]=t;
        time[size-1]=1;}
    void draw(GLenum type=GL_LINES){
        Vector p, c;
        bool allow=false;
        glBegin(type);
        for(int i=0;i<res;i++){
            double t=i/(res-1.);
            c=r(t);
            if(allow && canDraw(t) && size>3){
                Vector currentCopy=c;
                Vector prevCopy=p;
                if(cohenSutherland(prevCopy, currentCopy, lb, rt) && c.x==c.x){
                    glVertex2f(prevCopy.x, prevCopy.y);
                    glVertex2f(currentCopy.x, currentCopy.y); }}
            if(canDraw(t) && size>3) allow=true;
            p=c; }
        glEnd(); }
    bool pick(const Vector& lb, const Vector& rt){
        Vector prev, current;
        for(int i=0;i<res;i++){
            double t=i/(res-1.);
            if(canDraw(t)){
                current=r(t);
                Vector currentCopy=current;
                Vector prevCopy=prev;
                if(i>1 && cohenSutherland(prevCopy, currentCopy, lb, rt)) return true;
                prev=current; }}
        return false; }
    //Cohen-Sutherland vagasi algoritmus: http://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland_algorithm
    bool cohenSutherland(Vector& p0, Vector& p1, const Vector& lb, const Vector& rt){
        unsigned char c0=code(p0, lb, rt), c1=code(p1, lb, rt);
        Vector p;
        while(1){
            if(!(c0|c1)) return true;
            else if(c0&c1) return false;
            else {
                unsigned char c=c0?c0:c1;
                if(c&1){
                    p.x=p0.x+(p1.x-p0.x)*(rt.y-p0.y)/(p1.y-p0.y);
                    p.y=rt.y;
                } else if(c&2){
                    p.x=p0.x+(p1.x-p0.x)*(lb.y-p0.y)/(p1.y-p0.y);
                    p.y=lb.y;
                } else if(c&4){
                    p.y=p0.y+(p1.y-p0.y)*(rt.x-p0.x)/(p1.x-p0.x);
                    p.x=rt.x;
                } else {
                    p.y=p0.y+(p1.y-p0.y)*(lb.x-p0.x)/(p1.x-p0.x);
                    p.x=lb.x; }
                if(c==c0){
                    p0=p;
                    c0=code(p0, lb, rt);
                } else {
                    p1=p;
                    c1=code(p1, lb, rt); }}}}
    unsigned char code(Vector& p, const Vector& lb, const Vector& rt){ return (p.x<lb.x)<<3|(p.x>rt.x)<<2|(p.y<lb.y)<<1|(p.y>rt.y); }
    void addPoint(float x, float y){
        if(x>lb.x&&x<rt.x&&y>lb.y&&y<rt.y&&size<pMax){
            points[size++]=Vector(x,y);
            redistribute(); }}
    void move(Vector& v){ for(int i=0;i<size;i++){ points[i]=points[i]+v; }}
    virtual Vector r(double t)=0;
    virtual bool canDraw(double t)=0;
    protected:
    int getIndex(double t) { for(int i=0;i<size;i++) if(time[i+1]>t) return i; return -1; }};
class CRSpline: public Spline{
    public:
    CRSpline():Spline(){}
    Vector r(double t){
        int i=getIndex(t);
        t-=time[i];
        Vector* p=points;
        double* ti=time;
        double dt=time[i+1]-time[i];
        Vector v0=((p[i]-p[i-1])/(ti[i]-ti[i-1])+(p[i+1]-p[i])/(ti[i+1]-ti[i]))/2.;
        Vector v1=((p[i+1]-p[i])/(ti[i+1]-ti[i])+(p[i+2]-p[i+1])/(ti[i+2]-ti[i+1]))/2.;
        Vector a=(v1+v0)/dt/dt-(p[i+1]-p[i])*2/dt/dt/dt;
        Vector b=(p[i+1]-p[i])*3/dt/dt-(v1+v0*2)/dt;
        Vector c=v0;
        Vector d=p[i];
        return a*t*t*t+b*t*t+c*t+d; }
    bool canDraw(double t){ return getIndex(t)>0 && getIndex(t)<size-2; }};
class KKSpline: public Spline{
    public:
    KKSpline():Spline(){}
    Vector r(double t){
        int i=getIndex(t);
        double* ti=time;
        double tt=(t-ti[i])/(ti[i+1]-ti[i]);
        if(i==0) return lagrange(t, &(points[i]), 3, &(ti[i]));
        else if(i>size-3) return lagrange(t, &(points[size-3]), 3, &(ti[size-3]));
        Vector r1=lagrange(t, &(points[i-1]), 3, &(ti[i-1]));
        Vector r2=lagrange(t, &(points[i]), 3, &(ti[i]));
        return r1*(1-tt)+r2*(tt); }
    bool canDraw(double t)  {return getIndex(t)>-1; }
    protected:
    Vector lagrange(double t, Vector* p, int n, double* ti){
        Vector sum;
        for(int i=0;i<n;i++) {
            double prod=1;
            for(int j=0;j<n;j++) {
                if(i!=j){
                    prod*=(t-ti[j]);
                    prod/=(ti[i]-ti[j]); }}
            sum=sum+p[i]*prod; }
        return sum; }};
const int screenWidth = 600;
const int screenHeight = 600;
const double zoomf=10;
Vector dragStart;
Vector dragEnd;
bool dragging=false;
bool dragtype=false;
const int dzone=2;
Vector lb=Vector(100, 100);
Vector rt=Vector(500, 500);
KKSpline kk;
CRSpline cr;
void onInitialization( ) {
    glViewport(0, 0, screenWidth, screenHeight);
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(lb.x, rt.x, lb.y, rt.y); }
void onDisplay( ) {
    glClearColor(0, 0, 0, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3f(0,1,0);
    cr.draw();
    glColor3f(1,0,0);
    kk.draw();
    glutSwapBuffers(); }
void onKeyboard(unsigned char key, int x, int y) {
    double xx=((double)x/screenWidth);
    double yy=1-((double)y/screenHeight);
    xx=lb.x+xx*(rt.x-lb.x);
    yy=lb.y+yy*(rt.y-lb.y);
    if (key == 'z') {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        lb.x=xx+(lb.x-xx)/zoomf;
        rt.x=xx+(rt.x-xx)/zoomf;
        lb.y=yy+(lb.y-yy)/zoomf;
        rt.y=yy+(rt.y-yy)/zoomf;
        gluOrtho2D(lb.x, rt.x, lb.y, rt.y);
        glutPostRedisplay(); }
    if (key == 'Z') {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        lb.x=xx+(lb.x-xx)*zoomf;
        rt.x=xx+(rt.x-xx)*zoomf;
        lb.y=yy+(lb.y-yy)*zoomf;
        rt.y=yy+(rt.y-yy)*zoomf;
        gluOrtho2D(lb.x, rt.x, lb.y, rt.y);
        glutPostRedisplay(); }}
void onMouse(int button, int state, int x, int y) {
    double xx=((double)x/screenWidth);
    double yy=1-((double)y/screenHeight);
    xx=lb.x+xx*(rt.x-lb.x);
    yy=lb.y+yy*(rt.y-lb.y);
    if (button == GLUT_LEFT && state == GLUT_DOWN){
        kk.addPoint(xx, yy);
        cr.addPoint(xx, yy);
        glutPostRedisplay( ); }
    if(button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN && !dragging) {
        Vector mlb=Vector(lb.x+(x-dzone+0.)/screenWidth*(rt.x-lb.x), lb.y+(1-((y+dzone+0.)/screenHeight))*(rt.y-lb.y));
        Vector mrt=Vector(lb.x+(x+dzone+0.)/screenWidth*(rt.x-lb.x), lb.y+(1-((y-dzone+0.)/screenHeight))*(rt.y-lb.y));
        if(kk.pick(mlb, mrt)){
            dragging=true;
            dragtype=false;
            dragStart=Vector(xx,yy);
        }else if(cr.pick(mlb, mrt)){
            dragging=true;
            dragtype=true;
            dragStart=Vector(xx,yy);}}
    else if(button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN && dragging) {
        dragging=false;
        Vector end=Vector(xx,yy)-dragStart;
        if(dragtype) cr.move(end);
        else kk.move(end);
        glutPostRedisplay( ); }}
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
    return 0;}
