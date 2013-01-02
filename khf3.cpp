#include <math.h>
#include <stdlib.h>
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double min(double a, double b){
    return a<b?a:b;
}
double max(double a, double b){
    return a>b?a:b;
}
struct Vector {
   float x, y, z;
   Vector( ) {
    x = y = z = 0;
   }
   Vector(float x0, float y0, float z0 = 0) {
    x = x0; y = y0; z = z0;
   }
   Vector operator*(float a) {
    return Vector(x * a, y * a, z * a);
   }
   Vector operator/(float a) {
    return Vector(x / a, y / a, z / a);
   }
   Vector operator+(const Vector& v) {
    return Vector(x + v.x, y + v.y, z + v.z);
   }
   Vector operator-(const Vector& v) {
    return Vector(x - v.x, y - v.y, z - v.z);
   }
   float operator*(const Vector& v) {
    return (x * v.x + y * v.y + z * v.z);
   }
   Vector operator%(const Vector& v) {
    return Vector(y*v.z-z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
   }
   float length() { return sqrt(x * x + y * y + z * z); }
   float lengthSq() { return x*x+y*y+z*z; }
   Vector unit(){
       return (*this)*(1/length());
   }
};
struct Color {
   float r, g, b;
 
   Color( ) {
    r = g = b = 0;
   }
   Color(float r0, float g0, float b0) {
    r = r0; g = g0; b = b0;
   }
   Color operator*(float a) const{
    return Color(r * a, g * a, b * a);
   }
   Color operator*(const Color& c) const{
    return Color(r * c.r, g * c.g, b * c.b);
   }
   Color operator+(const Color& c) const{
    return Color(r + c.r, g + c.g, b + c.b);
   }
};
struct Photon:public Vector{
    Color c;
    Photon():Vector(),c(0,0,0){}
    Photon(Vector v, Color c):Vector(v),c(c){
        this->c.r=max(-c.r,c.r);
        this->c.g=max(-c.g,c.g);
        this->c.b=max(-c.b,c.b);
    }
};
struct Ray;
struct Obj;
struct Mat;
const int screenWidth = 600;
const int screenHeight = 600;
const int phms=150;
const int objLim=1000;
const long long phLim=10000;
const double EPS=1e-4;
const Color La(0.53,0.81,0.92);
Color lint(1,1,1);
Color image[screenWidth*screenHeight];
Obj* obj[objLim];
Color photonMap[phms*phms];
int objSize=0;
long long phSize=0;
Vector light(1.5,0,0.6);
struct Ray{
    Vector start;
    Vector dir;
    Ray(Vector& start, Vector& lookat):start(start),dir((lookat-start).unit()){}
    Vector operator[](double t){
        return start+dir*t;
    }
};
struct Mat{
    Color ka,kd,ks,n,k;
    double shine;
    bool reflective,refractive;
};
Mat dummyMat;
struct Obj{
    Mat mat;
    Obj(Mat& mat):mat(mat){};
    virtual double intersect(Ray& ray, Vector& normal)=0;
};
struct Sphere:public Obj{
    Vector center;
    double R;
    Sphere(Mat& mat,const Vector& center, double r):Obj(mat),center(center),R(r){}
    double intersect(Ray& ray, Vector& normal){
        double a=ray.dir*ray.dir*2;
        double b=ray.dir*(ray.start-center)*2;
        double c=(ray.start-center)*(ray.start-center)-R*R;
        double d=b*b-2*a*c;
        if(d<0) return -1;
        d=sqrt(d);
        double t1=(-b-d)/a;
        double t2=(-b+d)/a;
        if(t1<EPS&&t2<EPS) return -1;
        else if(t1<EPS||t2<EPS) {
            normal=(ray[max(t1,t2)]-center)/R;
            return max(t1,t2);
        }
        normal=(ray[min(t1,t2)]-center)/R;
        return min(t1,t2);
    }
};
struct Triangle:public Obj{
    Vector p1,p2,p3,n1,n2,n3,n;
    Triangle():Obj(dummyMat){};
    Triangle(Mat& mat, Vector& p1, Vector& p2,  Vector& p3,  Vector& n1,  Vector& n2,  Vector& n3):
    Obj(mat),p1(p1),p2(p2),p3(p3),n1(n1),n2(n2),n3(n3),n(((p2-p1)%(p3-p1)).unit()){};
    double intersect(Ray& ray, Vector& normal){
        double t=(p1-ray.start)*n/(ray.dir*n);
        if(t<EPS) return -1;
        Vector x=ray[t];
        if((p2-p1)%(x-p1)*n>-EPS&&
           (p3-p2)%(x-p2)*n>-EPS&&
           (p1-p3)%(x-p3)*n>-EPS){
            normal=n;
            return t;
           }
        return -1;
    }
};
struct Torus:public Obj{
    Triangle triangles[98];
    Sphere boundVol;
    double R,r;
    Torus(Mat& mat, double R, double r):Obj(mat),R(R),r(r),boundVol(dummyMat,Vector(0,0,0),0){
        boundVol=Sphere(mat, Vector(0,0,0), R+r);
        int us=7;
        int vs=7;
        double du=3.1415*2/7;
        double dv=3.1415*2/7;
        int nn=0;
        for(int u=0;u<us;u++){
            for(int v=0;v<vs;v++){
                Vector p00=p(u*du, v*dv);
                Vector p10=p((u+1)*du, v*dv);
                Vector p01=p(u*du, (v+1)*dv);
                Vector p11=p((u+1)*du, (v+1)*dv);
                Vector n00=n(u*du, v*dv);
                Vector n10=n((u+1)*du, v*dv);
                Vector n01=n(u*du, (v+1)*dv);
                Vector n11=n((u+1)*du, (v+1)*dv);
                triangles[nn++]=Triangle(mat, p00, p01, p10, n00, n01, n10);
                triangles[nn++]=Triangle(mat, p11, p10, p01, n11, n10, n01);
            }
        }
    }
    double intersect(Ray& r, Vector& normal){
        if(boundVol.intersect(r, normal)<0) return -1;
        double t=1e30;
        Vector norm;
        bool found=false;
        int in=-1;
        for(int i=0;i<98;i++){
            double t2=triangles[i].intersect(r, norm);
            if(t2<t&&t2>EPS){
                found=true;
                t=t2;
                in=i;
            }
        }
        if(found){
            Vector a=triangles[in].p3-triangles[in].p1;
            Vector b=triangles[in].p2-triangles[in].p1;
            Vector c=r[t]-triangles[in].p1;
            double aa=a*a;
            double ab=a*b;
            double ac=a*c;
            double bb=b*b;
            double bc=b*c;
            double det=1/(aa*bb-ab*ab);
            double u=(bb*ac-ab*bc)*det;
            double v=(aa*bc-ab*ac)*det;
            normal=triangles[in].n1*(1-u-v)+triangles[in].n2*(v)+triangles[in].n3*(u);
            return t;
        }
        return -1;
    }
    private:
    Vector p(double u, double v){
        return Vector((R+r*cos(v))*cos(u), r*sin(v),(R+r*cos(v))*sin(u));
    }
    Vector n(double u, double v){
        return Vector(cos(v)*cos(u), sin(v), cos(v)*sin(u)).unit();
    }
};
Vector reflectDir(Vector& v, Vector& n){
    return (v-n*(v*n)*2).unit();
}
Vector refractDir(Vector& v, Vector n, double ior, bool& valid){
    double tst=v*n;
    if(tst<0){
        n=n*-1;
        ior=1/ior;
    }
    double cost=v*n*1;
    double d=1-(1-cost*cost)*ior*ior;
    if(d<0){
        valid=false;
        return Vector(0,0,0);
    }
    return (v*ior-n*(cost*ior-sqrt(d))).unit();
}
double calcFresnel(Vector& v, Vector n, double ior, double k){
    double tst=v*n;
    if(tst<0){
        n=n*-1;
        ior=1/ior;
    }
    ior=1/ior;
    double f0=((ior-1)*(ior-1)+k*k)/((ior+1)*(ior+1)+k*k);
    return f0+(1-f0)*pow(1-n*v, 5);
}
Color fresnel(Vector& v, Vector &n, Mat& m, bool inv){
    double r=calcFresnel(v,n,m.n.r,m.k.r);
    double g=calcFresnel(v,n,m.n.g,m.k.g);
    double b=calcFresnel(v,n,m.n.b,m.k.b);
    if(inv)
        return Color(1-r, 1-g, 1-b);
    return Color(r,g,b);
}
double firstIntersect(Ray& ray, Vector& x, Vector& n, Obj** o){
    double t=1e30;
    Vector norm;
    bool found=false;
    for(int i=0;i<objSize;i++){
        double t2=obj[i]->intersect(ray, norm);
        if(t2<t&&t2>EPS){
            *o=obj[i];
            n=norm;
            found=true;
            t=t2;
        }
    }
    if(found){
        x=ray[t];
        return t;
    }
    return -1;
}
Color directLight(Vector& v, Vector& x, Vector& n, Mat& m){
    Color c(m.ka);
    Ray r(x, light);
    Vector t1,t2;
    Obj* t3;
    double dst=firstIntersect(r, t1, t2, &t3);
    if(dst<EPS||dst>(x-light).length()){
        c=c+lint*m.kd*max(0,n*(light-x).unit());
        c=c+lint*m.ks*max(0,pow(n*(light-x-v).unit(), m.shine));
    }
    return c;
}
Color& ph(int x, int y){
    return photonMap[y*phms+x];
}
Color ph2(int x, int y){
    int r=phms/25;
    if(x>=phms-r) x=phms-r;
    if(x<r) x=r;
    if(y>=phms-r) y=phms-r;
    if(y<r) y=r;
    Color sum;
    for(int i=-r;i<r+1;i++){
        for(int j=-r;j<r+1;j++){
            sum=sum+ph(x+i,y+j)*max(0,r-sqrt(i*i+j*j));
        }
    }
    return sum*(1./r/r);
}
Color trace(Ray& ray, int i){
    if(!i) return Color(0,0,0);
    Obj* o;
    Vector x;
    Vector n;
    if(firstIntersect(ray, x, n, &o)<0) return La;
    Color c=directLight(ray.dir,x,n,o->mat);
    if(!o->mat.reflective&&!o->mat.refractive){
        double px=(x.x/4+0.5)*phms;
        double py=(x.y/4+0.5)*phms;
        double u=px-(int)px;
        double v=py-(int)py;
        Color f=ph2((int)px, (int)py)*(1-u)*(1-v)+
            ph2((int)px+1, (int)py)*(u)*(1-v)+
            ph2((int)px, (int)py+1)*(1-u)*(v)+
            ph2((int)px+1, (int)py+1)*(u)*(v);
            c=c+f*(5000./phLim);
    }
    if(o->mat.reflective) {
        Vector refd=reflectDir(ray.dir,n)+x;
        Ray ref(x,refd);
        if(o->mat.refractive){
            bool valid=true;
            Vector refrd=refractDir(ray.dir,n,o->mat.n.r,valid);
            if(valid){
                c=c+fresnel(ray.dir,n,o->mat,false)*trace(ref,i-1);
            } else {
                c=c+trace(ref,i-1);
            }
        } else {
            c=c+fresnel(ray.dir,n,o->mat,false)*trace(ref,i-1);
        }
    }
    if(o->mat.refractive){
        bool valid=true;
        Vector refrd=refractDir(ray.dir,n,o->mat.n.r,valid)+x;
        if(valid){
            Ray refr(x,refrd);
            c=c+fresnel(ray.dir,n,o->mat,true)*trace(refr,i-1);
        }
    }
    return c;
}
void render(){
    Vector eye(0,2,1);
    Vector up(0,0,1);
    Vector lookat(0,0,0);
    double fov=80/180.*3.1415;
    Vector dir=(lookat-eye).unit();
    up=(up-dir*(dir*up)).unit();
    Vector right=up%dir;
    lookat=eye+dir/tan(fov/2);
    for(int y=0;y<screenHeight;y++){
        for(int x=0;x<screenWidth;x++){
            Vector look=lookat+up*(2.*y/screenHeight-1)+right*(2.*x/screenWidth-1);
            Ray ray(eye, look);
            image[y*screenWidth+x]=trace(ray, 8);
        }
    }
}
void add(Obj* o){
    if(objSize<objLim){
        obj[objSize++]=o;
    }
}
void addPhoton(const Vector& x, const Color& c){
    int px=(x.x/4+0.5)*phms;
    int py=(x.y/4+0.5)*phms;
    photonMap[py*phms+px]=photonMap[py*phms+px]+c;
    phSize++;
}
void shoot(const Color& pow, Ray& ray, int i){
    if(i>8) return;
    Obj* o;
    Vector x;
    Vector n;
    if(firstIntersect(ray, x, n, &o)<0) return;
    if(i&&!o->mat.reflective&&!o->mat.refractive){
        if(n*ray.dir<0){
            addPhoton(x, pow*(n*ray.dir*-1));
        }
        return;
    }
    if(o->mat.reflective) {
        Vector refd=reflectDir(ray.dir,n)+x;;
        Ray ref(x,refd);
        if(o->mat.refractive){
            bool valid=true;
            Vector refrd=refractDir(ray.dir,n,o->mat.n.r,valid);
            if(valid){
                shoot(pow*fresnel(ray.dir,n,o->mat,false), ref, i+1);
            } else {
                shoot(pow, ref, i+1);
            }
        } else {
            shoot(pow*fresnel(ray.dir,n,o->mat,false), ref, i+1);
        }
    }
    if(o->mat.refractive){
        bool valid=true;
        Vector refrd=refractDir(ray.dir,n,o->mat.n.r,valid)+x;
        if(valid){
            Ray refr(x,refrd);
            shoot(pow*fresnel(ray.dir,n,o->mat,true), refr, i+1);
        }
    }
}
void gen(long long n){
    for(long long i=0;i<n;i++){
        Vector v;
        do{
            v.x=(double)rand()/RAND_MAX-1;
            v.y=(double)rand()/RAND_MAX-0.5;
            v.z=(double)rand()/RAND_MAX-1;
        }while(v.length()>0.25&&fabs(v.y)>fabs(v.x));
        v=v.unit()+light;
        Ray r(light, v);
        shoot(Color(1,1,1), r, 0);
    }
}
double lum(Color& c){
    return 0.21*c.r+0.72*c.g+0.07*c.b;
}
void tone(){
    double maxl=0;
    double a=0.8;
    for(int y=0;y<screenHeight;y++){
        for(int x=0;x<screenWidth;x++){
            double l=lum(image[y*screenWidth+x]);
            maxl=maxl>l?maxl:l;
        }
    }
    for(int y=0;y<screenHeight;y++){
        for(int x=0;x<screenWidth;x++){
           double l=lum(image[y*screenWidth+x]);
           double d=l*a/(l*a-l+maxl);
           image[y*screenWidth+x]=image[y*screenWidth+x]*d;
        }
    }
    double t=1;
}
void onInitialization( ) {
    glViewport(0, 0, screenWidth, screenHeight);
    srand(42);
    Mat gold;
    gold.ka=Color(0.24725, 0.1995, 0.0745);
    gold.kd=Color(0.75164, 0.60648, 0.22648);
    gold.ks=Color(0.628281, 0.555802, 0.366065);
    gold.shine=0.4 * 128;
    gold.n=Color(0.17, 0.35, 1.5);
    gold.k=Color(3.1, 2.7, 1.9);
    gold.reflective=true;
    gold.refractive=false;
    Mat brown;
    brown.ka=Color(0.1,0.1,0.1);
    brown.kd=Color(0.58, 0.294, 0);
    brown.ks=Color(0,0,0);
    brown.shine=0;
    brown.reflective=false;
    brown.refractive=false;
    double gn=1.5;
    Mat glass;
    glass.ka=Color(0.1,0.1,0.1);
    glass.kd=Color(0.2,0.2,0.2);
    glass.ks=Color(1,1,1);
    glass.shine=120;
    glass.reflective=true;
    glass.refractive=true;
    glass.n=Color(gn,gn,gn);
    glass.k=Color(0,0,0);
    Mat glass2=glass;
    glass2.n=Color(1/gn, 1/gn, 1/gn);
    add(new Sphere(glass2, Vector(0,0,0), 0.4));
    add(new Torus(gold, 0.25, 0.1));
    double z=-0.65;
    Vector p00(-2,-2,z);
    Vector p10(2,-2,z);
    Vector p01(-2,2,z);
    Vector p11(2,2,z);
    Vector n(0,0,1);
    add(new Triangle(brown, p00, p10, p01, n,n,n));
    add(new Triangle(brown, p11, p01, p10, n,n,n));
    int cube[12][3][3]={
            {
                {-1, -1, -1},
                {1, -1, -1},
                {1, -1, 1}
            },
            {
                {-1, -1, 1},
                {-1, -1, -1},
                {1, -1, 1}
            },
            {
                {-1, 1, -1},
                {1, 1, 1},
                {1, 1, -1}
            },
            {
                {-1, 1, 1},
                {1, 1, 1},
                {-1, 1, -1}
            },
            {
                {-1, -1, -1},
                {-1, 1, 1},
                {-1, 1, -1}
            },
            {
                {-1, 1, 1},
                {-1, -1, -1},
                {-1, -1, 1}
            },
            {
                {1, -1, -1},
                {1, 1, -1},
                {1, 1, 1}
            },
            {
                {1, 1, 1},
                {1, -1, 1},
                {1, -1, -1}
            },
            {
                {-1, 1, -1},
                {1, 1, -1},
                {-1, -1, -1}
            },
            {
                {1, -1, -1},
                {-1, -1, -1},
                {1, 1, -1}
            },
            {
                {-1, 1, 1},
                {-1, -1, 1},
                {1, 1, 1},
            },
            {
                {1, -1, 1},
                {1, 1, 1},
                {-1, -1, 1}
            }};
    double a=0.5,b=0.5,c=0.5;
    for(int i=0;i<12;i++){
        Vector p1(cube[i][0][0]*a, cube[i][0][1]*b, cube[i][0][2]*c);
        Vector p2(cube[i][1][0]*a, cube[i][1][1]*b, cube[i][1][2]*c);
        Vector p3(cube[i][2][0]*a, cube[i][2][1]*b, cube[i][2][2]*c);
        Vector n=((p2-p1)%(p3-p1)).unit();
        add(new Triangle(glass, p1, p2, p3, n, n, n));
    }
    gen(phLim);
    render();
    tone();
    for(int i=0;i<objSize;i++){
        delete obj[i];
    }
}
void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);
    glutSwapBuffers();
 
}
void onKeyboard(unsigned char key, int x, int y) {
}
void onMouse(int button, int state, int x, int y) {
}
void onIdle( ) {
}
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