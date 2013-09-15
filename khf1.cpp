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
    Vector operator+(const Vector &v) { return Vector(x + v.x, y + v.y, z + v.z); }
    Vector operator-(const Vector &v) { return Vector(x - v.x, y - v.y, z - v.z); }
    float operator*(const Vector &v) { return (x * v.x + y * v.y + z * v.z); }
    Vector operator%(const Vector &v) { return Vector(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
    float Length() { return sqrt(x * x + y * y + z * z); }};
struct Color {
    float r, g, b;
    Color( ) { r = g = b = 0; }
    Color(float r0, float g0, float b0) { r = r0; g = g0; b = b0; }
    Color operator*(float a) { return Color(r * a, g * a, b * a); }
    Color operator*(const Color &c) { return Color(r * c.r, g * c.g, b * c.b); }
    Color operator+(const Color &c) { return Color(r + c.r, g + c.g, b + c.b); }};
struct Actor;
struct Tower;
float f(float, float);
float gradf(float, float, Vector);
float col(float);
const int screenWidth = 600, screenHeight = 600;
long lastSim = 0, dt = 100, frame = 0, frameCanTransmit = 0;
Color image[screenWidth *screenHeight];
struct Actor {
    static const float vmax = 0.02;
    Color c;
    Vector r, v, d;
    bool isBoy;
    Actor(Color c, bool isBoy): c(c), isBoy(isBoy) {};
    void draw();
    void redirect(float x, float y);
    void tick();
    bool canTransmit(Tower t);};
struct Tower {
    Vector r;
    void relocate();
    void draw();};
Actor hansel(Color(0.16, 0.556, 1.0), true), gretel(Color(1.0, 0.67, 1.0),false);
Tower tower;
void Actor::draw() {
    float R=0.05, yOffset=0.05, xOffset=0;
    glColor3f(c.r, c.g, c.b);
    glBegin(GL_LINE_STRIP); for(float i=0;i<3.1415*2;i+=3.1415*2/100) glVertex2f(r.x+xOffset+R*sin(i), r.y+yOffset+R*cos(i)); glEnd();
    if(isBoy) {
        glBegin(GL_LINES);
        glVertex2f(r.x+xOffset+R*sin(3.1415/4), r.y+yOffset+R*cos(3.1415/4));
        glVertex2f(r.x+xOffset+R*2, r.y+yOffset+R*2);
        glVertex2f(r.x+xOffset+R*2, r.y+yOffset+R*2);
        glVertex2f(r.x+xOffset+R, r.y+yOffset+R*2);
        glVertex2f(r.x+xOffset+R*2, r.y+yOffset+R*2);
        glVertex2f(r.x+xOffset+R*2, r.y+yOffset+R);
        glEnd();
    } else {
        glBegin(GL_LINES);
        glVertex2f(r.x+xOffset, r.y+yOffset-R);
        glVertex2f(r.x+xOffset, r.y+yOffset-R*3);
        glVertex2f(r.x+xOffset-R, r.y+yOffset-R*2);
        glVertex2f(r.x+xOffset+R, r.y+yOffset-R*2);
        glEnd(); }}
void Actor::redirect(float nx, float ny) {
    float dx = r.x - nx, dy = r.y - ny, vx = -sin(atan2(dx, dy)) * vmax, vy = -cos(atan2(dx, dy)) * vmax;
    d.x = vx;
    d.y = vy;}
void Actor::tick() {
    if (d.Length() < 0.001) {
        r.z = f(r.x, r.y); return; }
    float m = atan(gradf(r.x, r.y, d) * 0.0002) / 3.14 * 180;
    v.x = d.x * (1 - m / 90.0) * cos(m * 3.14 / 180);
    v.y = d.y * (1 - m / 90.0) * cos(m * 3.14 / 180);
    r.x += v.x;
    r.y += v.y;
    r.z = f(r.x, r.y);
    if (r.x > 1 || r.x < -1) {
        d.x *= -1;
        if (r.x > 1) r.x = 1;
        if (r.x < -1) r.x = -1;
    } if (r.y > 1 || r.y < -1) {
        d.y *= -1;
        if (r.y > 1) r.y = 1;
        if (r.y < -1) r.y = -1; }}
bool Actor::canTransmit(Tower t) {
    Vector end = t.r + Vector(0, 0, 20), trace = r - end;
    for (float i = 0; i < 1; i += 0.001) {
        Vector tt = end + (trace * i);
        if (f(tt.x, tt.y) > tt.z) return false; }
    return true; }
void Tower::relocate() {
    r.x = (float)rand() / RAND_MAX * 2.0 - 1.0;
    r.y = (float)rand() / RAND_MAX * 2.0 - 1.0;
    r.z = f(r.x, r.y); }
void Tower::draw() {
    glColor3f(1.0, 1.0, 0.0);
    glBegin(GL_TRIANGLES);
    glVertex2f(r.x - 0.05, r.y - 0.05); glVertex2f(r.x, r.y + 0.05); glVertex2f(r.x + 0.05, r.y - 0.05);
    glEnd(); }
void onInitialization() {
    glViewport(0, 0, screenWidth, screenHeight);
    srand(glutGet(GLUT_ELAPSED_TIME));
    tower.relocate();
    for (int y = 0; y < screenHeight; y++) for (int x = 0; x < screenWidth; x++) {
            float z=col(f((float)x/screenWidth*2-1, (float)y/screenHeight*2-1));
            image[y*screenWidth+x]=Color(z,z,z); }}
float f(float x, float y) {
    return 250.0 + (cos(x * y * 10) * sin(x * 6) + 1) * (1014 - 250) / 2.0; }
float gradf(float x, float y, Vector v) {
    if(v.Length() < 0.0001) return 0;
    v = v * (1 / v.Length());
    return Vector( 2292 * cos(6 * x) * cos(10 * x * y) - 3820 * y * sin(6 * x) * sin(10 * x * y), -3820 * x * sin(6 * x) * sin(10 * x * y))*v; }
float col(float h) {
    return (h - 250.0) / (1014.0 - 250.0); }
void onDisplay() {
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);
    hansel.draw();
    gretel.draw();
    tower.draw();
    if (hansel.canTransmit(tower) && gretel.canTransmit(tower)) {
        glColor3f(1.0, 0.0, 0.0);
        glBegin(GL_TRIANGLES);
        glVertex2f(hansel.r.x, hansel.r.y); glVertex2f(gretel.r.x, gretel.r.y); glVertex2f(tower.r.x, tower.r.y);
        glEnd();
    }
    glColor3f(.0, .0, .0);
    glBegin(GL_POLYGON);
    glVertex2f(-.5, .95); glVertex2f(.5, .95); glVertex2f(.5, .85); glVertex2f(-.5, .85);
    glEnd();
    glColor3f(.0, 1., .0);
    glBegin(GL_POLYGON);
    glVertex2f(-.5, .95); glVertex2f((float)frameCanTransmit / frame - 0.5, .95); glVertex2f((float)frameCanTransmit / frame - 0.5, .85); glVertex2f(-.5, .85);
    glEnd();
    glutSwapBuffers(); }
void onKeyboard(unsigned char key, int x, int y) {
    if (key == 't') { tower.relocate(); glutPostRedisplay(); }}
void onMouse(int button, int state, int x, int y) {
    Actor *e;
    if (button == GLUT_LEFT && state == GLUT_DOWN) e = &hansel;
    else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) e = &gretel;
    else return;
    float xx = (float)(x - screenWidth / 2.0) / screenWidth * 2.0;
    float yy = (screenHeight / 2.0 - (float)y) / screenHeight * 2.0;
    e->redirect(xx, yy); }
void onIdle() {
    long time = glutGet(GLUT_ELAPSED_TIME), t;
    if ((time - lastSim) / dt) {
        for (t = 0; t < (time - lastSim) / dt; t++) {
            frame++;
            hansel.tick();
            gretel.tick();
            if (gretel.canTransmit(tower) && hansel.canTransmit(tower)) frameCanTransmit++; }
        lastSim += t * dt;
        glutPostRedisplay(); }}
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
