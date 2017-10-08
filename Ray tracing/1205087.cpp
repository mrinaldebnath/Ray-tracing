#include <bits/stdc++.h>
#include "1205087.hpp"

#ifdef __linux
#include <GL/glut.h>
#else
#include <windows.h>
#include <glut.h>
#endif // windows

#define pi (2*asin(1.0))

using namespace std;



const int INF = 1000000006;

const double EPS = 1e-11;


struct Vector
{
    double x,y,z,homo;
    Vector()
    {
        x=y=z=0.0;
        homo=1.0;
    }
    Vector(double a,double b,double c)
    {
        x=a;
        y=b;
        z=c;
        homo=1.0;
    }

    Vector crossVector(Vector v)
    {
        Vector ret;
        ret.x=y*v.z-z*v.y;
        ret.y=-(x*v.z-z*v.x);
        ret.z=x*v.y-y*v.x;
        ret.homo=1.0;

        return ret;
    }


    void normalise()
    {


        double val=sqrt(x*x+y*y+z*z);
        if(val==0.0)return;

        x=x*1.0/val;
        y=y*1.0/val;
        z=z*1.0/val;
        homo=1.0;

    }
    Vector scale(double mul)
    {
        Vector ret;
        ret.x=x;
        ret.y=y;
        ret.z=z;
        ret.x*=1.0*mul;
        ret.y*=1.0*mul;
        ret.z*=1.0*mul;
        ret.homo*=1.0;
        return ret;
    }

    Vector minusVector(Vector b)
    {
        Vector ret;
        ret.x=x-b.x;
        ret.y=y-b.y;
        ret.z=z-b.z;
        ret.homo*=1.0;
        return ret;
    }

    Vector plusVector(Vector b)
    {
        Vector ret;
        ret.x=x+b.x;
        ret.y=y+b.y;
        ret.z=z+b.z;
        ret.homo*=1.0;
        return ret;
    }

    double dotProduct(Vector b)
    {
        double ret=0.0;
        ret+=x*b.x;
        ret+=y*b.y;
        ret+=z*b.z;
        return ret;
    }

};

double area(Vector a, Vector b, Vector c)
{
    Vector ac=a.minusVector(c);
    Vector bc=b.minusVector(c);

    Vector ret=ac.crossVector(bc);
    return 0.5*(sqrt(ret.dotProduct(ret)));
}

bool inside_point(Vector m,Vector p,Vector q,Vector r)
{
    double mpq=area(m,p,q);
    double mqr=area(m,q,r);
    double mrp=area(m,r,p);

    double pqr=area(p,q,r);
    double tot=mpq+mqr+mrp;

    if(abs(tot-pqr)<EPS) return true;
    return false;

}

void specialKeyListener(int key, int x,int y);
void keyboardListener(unsigned char key, int x,int y);
void init();

void drawAxes(double l);
void drawSquare(double x,double y,double z,double a,double c);
void drawSphere(double radius,int slices,int stacks,double x,double y,double z);
void printImage();
void rayTracer();
void drawCube(Vector centre);


struct Ray
{
    Vector p0;
    Vector dir;
    Ray()  {}
    Ray(Vector p0, Vector p1)
    {
        Vector temp=p0.minusVector(p1);
        dir = temp.scale(-1);
        dir.normalise();
    }

};


struct Color
{
    double r,g,b;
    Color() {}
    Color(double red,double green,double blue)
    {
        r=red;
        g=green;
        b=blue;
    }

    Color scale(double mul)
    {
        Color ret;
        ret.r=r;
        ret.g=g;
        ret.b=b;
        ret.r*=1.0*mul;
        ret.g*=1.0*mul;
        ret.b*=1.0*mul;
        return ret;
    }


    Color plusColor(Color c)
    {
        Color ret;
        ret.r=r+c.r;
        ret.g=g+c.g;
        ret.b=b+c.b;
        return ret;
    }

    Color dotProduct(Color c)
    {
        Color ret;
        ret.r=r*c.r;
        ret.g=g*c.g;
        ret.b=b*c.b;
        return ret;
    }

};

double eff=1.0;

struct Triangle
{
    Vector x,y,z;
    Vector normal;
    double a1,b1,c1,d1;
    Color c;
    double amb, diffuse, spec, refl;
    double shine;
    double refract_coeff;

    Triangle()  {}
    Triangle(Vector p0,Vector p1,Vector p2)
    {
        x=p0;
        y=p1;
        z=p2;

        trianglePlane();
    }

    void trianglePlane()
    {
        Vector dir1(y.x-x.x,y.y-x.y,y.z-x.z);
        Vector dir2(z.x-x.x,z.y-x.y,z.z-x.z);
        Vector n=dir1.crossVector(dir2);

        a1=n.x;
        b1=n.y;
        c1=n.z;
        normal=Vector(n.x,n.y,n.z);
        d1=-(a1*y.x+b1*y.y+c1*y.z);
    }


    double rayIntersection(Ray ri, Vector &n)         //returns the parametric value of t
    {
        ri.dir.normalise();

        double denom=a1*ri.dir.x+b1*ri.dir.y+c1*ri.dir.z;
        if(abs(denom)<EPS)return -INF;
        double numerator=-(d1+a1*ri.p0.x+b1*ri.p0.y+c1*ri.p0.z);

        double t=numerator/denom;

        Vector temp=ri.dir.scale(t);

        bool onTrngl = inside_point(ri.p0.plusVector(temp),x,y,z);
        if(onTrngl==false)  return -INF;
        n =normal;
        n.normalise();
        return t;
    }

    void draw()
    {
        glColor3f(c.r,c.g,c.b);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(z.x,z.y,z.z);
            glVertex3f(x.x,x.y,x.z);
            glVertex3f(y.x,y.y,y.z);
        }
        glEnd();

    }

};

struct Sphere
{
    double amb, diffuse, spec, refl;
    double rad;
    double shine;
    Vector center;
    double A,B,C,D,E,F,G,H,I,J;

    double refract_coeff;
    Color c;

    Sphere()  {}


    Vector rayIntersection(Ray ri, double &t)
    {
        double a=(A*ri.dir.x*ri.dir.x)+(B*ri.dir.y*ri.dir.y)+(C*ri.dir.z*ri.dir.z)+(D*ri.dir.x*ri.dir.y)+(E*ri.dir.y*ri.dir.z)+(F*ri.dir.z*ri.dir.x);
        double b=(2*A*ri.p0.x*ri.dir.x)+(2*B*ri.p0.y*ri.dir.y)+(2*C*ri.p0.z*ri.dir.z)+((D*ri.p0.x*ri.dir.y)+(D*ri.dir.x*ri.p0.y))+((E*ri.p0.y*ri.dir.z)+(E*ri.dir.y*ri.p0.z))+((F*ri.p0.x*ri.dir.z)+(F*ri.dir.x*ri.p0.z))+(G*ri.dir.x)+(H*ri.dir.y)+(I*ri.dir.z);
        double c=(A*ri.p0.x*ri.p0.x)+(B*ri.p0.y*ri.p0.y)+(C*ri.p0.z*ri.p0.z)+(D*ri.p0.x*ri.p0.y)+(E*ri.p0.y*ri.p0.z)+(F*ri.p0.z*ri.p0.x)+(G*ri.p0.x)+(H*ri.p0.y)+(I*ri.p0.z)+J;
        double d=(b*b)-(4*a*c);

        if(d<0)
        {
            t=-INF;
            return Vector(1,1,1);
        }

        d = sqrt(d);
        double t1 = (-b+d)/(2*a);
        double t2 = (-b-d)/(2*a);

        if(t1<0 && t2<0)
        {
            t=-INF;
            return Vector(1,1,1);
        }

        if(t1>=0 && t2>=0)
        {
            t=min(t1,t2);
        }
        else if(t1<0)t=t2;
        else t=t1;

        Vector tempS=ri.dir.scale(t);

        Vector temp2=ri.p0.plusVector(tempS);
        Vector normal = temp2.minusVector(center);
        normal.normalise();
        return normal;
    }

    void draw()
    {
        glPushMatrix();
        {
            glColor3f(c.r,c.g,c.b);
            drawSphere(rad, 60, 60,center.x, center.y,center.z);
        }
        glPopMatrix();

    }

    void read()
    {
        double xx,yy,zz;
        cin>>xx>>yy>>zz;
        cin>>rad;
        center=Vector(xx,yy,zz);

        cin>>xx>>yy>>zz;
        c=Color(xx,yy,zz);
        scanf("%lf %lf %lf %lf %lf",&amb,&diffuse,&spec,&refl,&shine);
        refract_coeff=min(refl*eff,1.0);

        A=1.0;
        B=1.0;
        C=1.0;
        D=0.0;
        E=0.0;
        F=0.0;
        G=-(2*center.x);
        H=-(2*center.y);
        I=-(2*center.z);
        J=center.dotProduct(center);
        J=J-(rad*rad);

    }
};


struct General
{
    Vector cube_reference;
    double A,B,C,D,E,F,G,H,I,J;
    double l,w,h;
    double refract_coeff;
    Color c;

    double amb, diffuse, spec, refl;
    double shine;

    General()   {}
    General(double a,double b,double c,double d,double e,double f,double g,double h,double i,double j)
    {
        A=a;
        B=b;
        C=c;
        D=d;
        E=e;
        F=f;
        G=g;
        H=h;
        I=i;
        J=j;
    }

    Vector get_normal(Vector intersection_point)
    {
        Vector ret;
        ret.x=(2*A*intersection_point.x)+(D*intersection_point.y)+(F*intersection_point.z)+G;
        ret.y=(2*B*intersection_point.y)+(D*intersection_point.x)+(E*intersection_point.z)+H;
        ret.z=(2*C*intersection_point.z)+(E*intersection_point.y)+(F*intersection_point.x)+I;

        ret.normalise();
        return ret;

    }

    Vector rayIntersection(Ray rri, double &t)
    {
        t = INF;
        Vector normal;

        double a=(A*rri.dir.x*rri.dir.x)+(B*rri.dir.y*rri.dir.y)+(C*rri.dir.z*rri.dir.z)+(D*rri.dir.x*rri.dir.y)+(E*rri.dir.y*rri.dir.z)+(F*rri.dir.z*rri.dir.x);
        double b=(2*A*rri.p0.x*rri.dir.x)+(2*B*rri.p0.y*rri.dir.y)+(2*C*rri.p0.z*rri.dir.z)+((D*rri.p0.x*rri.dir.y)+(D*rri.dir.x*rri.p0.y))+((E*rri.p0.y*rri.dir.z)+(E*rri.dir.y*rri.p0.z))+((F*rri.p0.x*rri.dir.z)+(F*rri.dir.x*rri.p0.z))+(G*rri.dir.x)+(H*rri.dir.y)+(I*rri.dir.z);
        double c=(A*rri.p0.x*rri.p0.x)+(B*rri.p0.y*rri.p0.y)+(C*rri.p0.z*rri.p0.z)+(D*rri.p0.x*rri.p0.y)+(E*rri.p0.y*rri.p0.z)+(F*rri.p0.z*rri.p0.x)+(G*rri.p0.x)+(H*rri.p0.y)+(I*rri.p0.z)+J;
        double d=(b*b)-(4*a*c);

        if(d<0)
        {
            t=-INF;
            normal=Vector(1,1,1);
            return normal;
        }
        d = sqrt(d);
        double t1 = (-b+d)/(2*a);
        double t2 = (-b-d)/(2*a);



        if(t1<0 && t2<0)
        {
            t=-INF;
            normal=Vector(1,1,1);
            return normal;
        }

        else if(t1<0 || t2<0)
        {
            t=max(t1,t2);

            Vector tempS=rri.dir.scale(t);

            Vector incident_pt=rri.p0.plusVector(tempS);

            bool inside_x=false;
            bool inside_y=false;
            bool inside_z=false;

            if(cube_reference.x<=incident_pt.x && (cube_reference.x+l)>=incident_pt.x)inside_x=true;
            if(cube_reference.y<=incident_pt.y && (cube_reference.y+w)>=incident_pt.y)inside_y=true;
            if(cube_reference.z<=incident_pt.z && (cube_reference.z+h)>=incident_pt.z)inside_z=true;

            if(l==0.0)inside_x=true;
            if(w==0.0)inside_y=true;
            if(h==0.0)inside_z=true;

            if(inside_x && inside_y && inside_z)

            {
                //cout<<"yes"<<endl;

                normal=get_normal(incident_pt);

                return normal;
            }

            else
            {
                normal=Vector(1,1,1);
                t=-INF;
                return normal;
            }

        }

        else
        {
            Vector tempS=rri.dir.scale(t1);

            Vector incident_pt=rri.p0.plusVector(tempS);

            bool inside_x=false;
            bool inside_y=false;
            bool inside_z=false;

            if(cube_reference.x<=incident_pt.x && (cube_reference.x+l)>=incident_pt.x)inside_x=true;
            if(cube_reference.y<=incident_pt.y && (cube_reference.y+w)>=incident_pt.y)inside_y=true;
            if(cube_reference.z<=incident_pt.z && (cube_reference.z+h)>=incident_pt.z)inside_z=true;

            if(l==0.0)inside_x=true;
            if(w==0.0)inside_y=true;
            if(h==0.0)inside_z=true;

            if(!inside_x || !inside_y || !inside_z)
            {
                t1=INF;
            }

            tempS=rri.dir.scale(t2);

            incident_pt=rri.p0.plusVector(tempS);

            inside_x=false;
            inside_y=false;
            inside_z=false;

            if(cube_reference.x<=incident_pt.x && (cube_reference.x+l)>=incident_pt.x)inside_x=true;
            if(cube_reference.y<=incident_pt.y && (cube_reference.y+w)>=incident_pt.y)inside_y=true;
            if(cube_reference.z<=incident_pt.z && (cube_reference.z+h)>=incident_pt.z)inside_z=true;

            if(l==0.0)inside_x=true;
            if(w==0.0)inside_y=true;
            if(h==0.0)inside_z=true;

            if(!inside_x || !inside_y || !inside_z)
            {
                t2=INF;
            }
            if(t1==INF && t2==INF)
            {
                t=-INF;
                normal=Vector(1,1,1);
                return normal;
            }
            else
            {
                t=min(t1,t2);

                Vector tempS=rri.dir.scale(t);

                Vector incident_pt=rri.p0.plusVector(tempS);

                normal=get_normal(incident_pt);

                return normal;
            }
        }


    }

};
Color **checkImage;
Vector cameraPos,cameraUp,cameraLook,cameraRight;
double screenHeight, screenWidth;
double recursionLevel;
bool texture;
Color **imageMap;
Color sourcePower;

struct CheckBoard
{
    double height, width;
    double boardHeight, boardWidth;
    double amb, diffuse, spec, refl;
    double shine;
    double refract_coeff;

    CheckBoard()  {}
    CheckBoard(double h, double w)
    {
        height = h;
        width = w;
        boardHeight = 25;
        boardWidth = 25;

        amb = .15;
        spec = .4,
        diffuse=.2,
        refl=.25;
        shine = 4;
        refract_coeff=.5;

    }


    bool onBoard(Vector p)
    {

        if(abs(p.x)<=width/2 && abs(p.y)<=height/2 && abs(p.z)<EPS) return true;
        return false;
    }

    Color getColor(Vector p)
    {
        if(!onBoard(p)) return Color(0,0,0);

        int row = (p.x+width/2), colm = (p.y+height/2);


        return checkImage[row][colm];

    }
    void draw()
    {
        double c;
        int colm = 0;
        for(double xx = (-width/2)+(boardWidth/2); xx< width/2; xx+= boardWidth, colm++)
        {
            int row = 0;
            for(double yy = (-height/2)+(boardHeight/2); yy<height/2; yy+= boardHeight, row++)
            {
                c=0.0;
                if((row+colm)%2)    c=1.0;
                drawSquare(xx,yy,0,boardWidth, c);
            }
        }
    }

    void draw2()
    {
        drawSquare(0,0,0,screenHeight,1);
    }

    Vector rayIntersection(Ray ri, double &t)
    {
        if(ri.dir.z<EPS && ri.dir.z>-EPS)
        {
            t=-INF;
            return Vector(1,1,1);
        }

        t = ri.p0.z/ri.dir.z;
        Vector temp=ri.dir.scale(-t);
        Vector pt = ri.p0.plusVector(temp);
        bool valid = onBoard(pt);
        t=-t;

        if(valid)
        {
            return Vector(0,0,1);
        }

        t=-INF;
        return Vector(1,1,1);

    }
};




void changeBoard()
{
    bitmap_image image1=bitmap_image("input.bmp");

    int w=image1.width();
    int h=image1.height();


    for(int i=0; i<(int)floor(screenWidth); i++)
    {
        for(int j=0; j<(int)floor(screenHeight); j++)
        {

            unsigned char r1,g1,b1;

            double i_double=(i*w)/screenWidth;
            double j_double=(j*h)/screenHeight;

            int i_=floor(i_double);
            int j_=floor(j_double);


            image1.get_pixel(i_,j_,r1,g1,b1);

            double r,g,b;
            r=r1/255.0;
            g=g1/255.0;
            b=b1/255.0;

            checkImage[i][(int)floor(screenHeight)-1-j]=Color(r,g,b);
        }
    }

}


Ray calculateReflection(Ray ri, Vector normal, Vector p0)
{
    Vector nm=normal.scale(-1);
    Vector new_ri=ri.dir.scale(-1);

    nm.normalise();
    new_ri.normalise();

    Ray ret;

    double t=new_ri.dotProduct(nm);
    double tempD=2*t;
    Vector temp=nm.scale(tempD);

    ret.dir = temp.minusVector(new_ri);
    ret.dir.normalise();
    ret.p0 = p0;
    return ret;
}


Ray calculateRefraction(Ray ri, Vector normal, Vector p0)
{


    Ray ret;
    ret.p0 = p0;
    ret.dir=ri.dir;
    ret.dir.normalise();
    return ret;
}

vector<Sphere>spheres;
vector<General>generals;
vector<Triangle>triangles;
vector<Vector>lights;
CheckBoard checkBoard;



#define CHECKBOARD 0
#define SPHERE 1
#define TRIANGLE 2
#define GENERAL 3
#define EYE 4


bool hasObstacle(Vector light, Vector pt)
{

    double tt=sqrt((light.x-pt.x)*(light.x-pt.x)+(light.y-pt.y)*(light.y-pt.y)+(light.z-pt.z)*(light.z-pt.z));
    Ray lightRay(light,pt);

    Vector tmp;
    double t;

    int s=spheres.size();
    for(int i = 0; i < s; i++)
    {
        tmp = spheres[i].rayIntersection(lightRay, t);

        if(t>EPS && tt>t+EPS)return true;

    }


    s=generals.size();
    for(int i = 0; i < s; i++)
    {
        tmp = generals[i].rayIntersection(lightRay, t);

        if(t>EPS && tt>t+EPS)return true;

    }

    s=triangles.size();
    for(int i = 0; i < s; i++)
    {
        t = triangles[i].rayIntersection(lightRay, tmp);

        if(t>EPS && tt>t+EPS)return true;

    }

    tmp = checkBoard.rayIntersection(lightRay, t);

    if(t>EPS && tt>t+EPS)return true;

    return false;
}

struct ObjectID
{
    int id;
    int type;
    ObjectID() {}
    ObjectID(int i)
    {
        type=EYE;
        id=i;
    }
    ObjectID(int i,int typ)
    {
        type=typ;
        id=i;
    }
    bool equals(ObjectID o)
    {
        return type==o.type && id==o.id;
    }
};


Color rayCast(Ray ri, int level, ObjectID from)
{
    if(level==recursionLevel)
    {
        return Color(0,0,0);
    }

    double spPt = INF, trPt= INF;
    double gnPt=INF;
    int sphr = -1, trnl=-1;
    int gnrl=-1;
    Vector spNormal,trNormal;
    Vector gnNormal;
    double amb, spec, diffuse, reflCoeff,shine;
    double refract_coeff;

    double tmp;
    Vector nm;

    int s=spheres.size();
    for(int i = 0; i < s; i++)
    {
        ObjectID ob=ObjectID(i,SPHERE);
        if(from.equals(ob)) continue;
        nm = spheres[i].rayIntersection(ri,tmp);


        if(tmp>0 && tmp<spPt)
        {
            sphr = i;
            spPt = tmp;
            spNormal = nm;
        }

    }
    if(spPt==INF)    sphr= -1;



    s=generals.size();
    for(int i = 0; i < s; i++)
    {
        ObjectID ob=ObjectID(i,GENERAL);
        if(from.equals(ob)) continue;
        nm = generals[i].rayIntersection(ri,tmp);

        if(tmp>0 && tmp<gnPt)
        {
            //cout<<"yes"<<endl;
            gnrl = i;
            gnPt = tmp;
            gnNormal = nm;
        }
    }
    if(gnPt==INF)    gnrl = -1;



    s=triangles.size();
    for(int i = 0; i < s; i++)
    {
        ObjectID ob=ObjectID(i,TRIANGLE);
        if(from.equals(ob)) continue;
        tmp = triangles[i].rayIntersection(ri,nm);
        if(tmp>0 && tmp<trPt)
        {
            trnl = i;
            trPt = tmp;
            trNormal = nm;
        }

    }
    if(trPt==INF)    trnl = -1;

    double chkPt;
    nm = checkBoard.rayIntersection(ri,chkPt);
    if(from.type==CHECKBOARD)
    {
        chkPt = -INF;
    }

    int selected;
    ObjectID newObj;

    Vector normal;
    Color colr;
    double minDist;

    if(sphr<0 && trnl<0 && gnrl<0)        //no intersection
    {
        if(chkPt<0 || chkPt>INF-EPS)
        {
            return Color(0,0,0);
        }
        selected = CHECKBOARD;
        normal = Vector(0,0,1);
        minDist = chkPt;
    }

    else if(trPt<=spPt && trPt<=gnPt) //triangle is closer
    {
        selected = TRIANGLE;
        minDist = trPt;
        if(chkPt>=0 && chkPt<trPt)
        {
            selected = CHECKBOARD;
            normal = Vector(0,0,1);
            minDist = chkPt;
        }
    }
    else if(spPt<=trPt && spPt<=gnPt)  //sphere is closer
    {
        selected = SPHERE;
        minDist = spPt;
        if(chkPt>=0 && chkPt<spPt)
        {
            selected = CHECKBOARD;
            normal = Vector(0,0,1);
            minDist = chkPt;
        }
    }

    else if(gnPt<=trPt && gnPt<=spPt)  //general is closer
    {
        //cout<<"yes"<<endl;
        selected = GENERAL;
        minDist = gnPt;
        if(chkPt>=0 && chkPt<gnPt)
        {
            selected = CHECKBOARD;
            normal = Vector(0,0,1);
            minDist = chkPt;
        }
    }


    Vector temp=ri.dir.scale(minDist);
    Vector incidentVector = ri.p0.plusVector(temp);

    if(selected == CHECKBOARD)
    {
        colr = checkBoard.getColor(incidentVector);
        amb = checkBoard.amb;
        spec = checkBoard.spec;
        diffuse=checkBoard.diffuse;
        reflCoeff=checkBoard.refl;
        shine = checkBoard.shine;
        refract_coeff=checkBoard.refract_coeff;
        newObj = ObjectID(0,CHECKBOARD);
    }
    else if(selected == GENERAL)
    {
        //cout<<"yes"<<endl;
        colr = generals[gnrl].c;
        amb = generals[gnrl].amb;
        spec = generals[gnrl].spec;
        diffuse = generals[gnrl].diffuse;
        reflCoeff = generals[gnrl].refl;
        shine = generals[gnrl].shine;
        normal = gnNormal;
        refract_coeff=generals[gnrl].refract_coeff;
        newObj = ObjectID(gnrl,GENERAL);
    }

    else if(selected == TRIANGLE)
    {
        diffuse = triangles[trnl].diffuse;
        amb = triangles[trnl].amb;
        colr = triangles[trnl].c;
        shine = triangles[trnl].shine;
        spec = triangles[trnl].spec;
        reflCoeff = triangles[trnl].refl;
        refract_coeff=triangles[trnl].refract_coeff;
        normal = trNormal;
        newObj = ObjectID(trnl,TRIANGLE);
    }

    else if(selected == SPHERE)
    {
        colr = spheres[sphr].c;
        amb = spheres[sphr].amb;
        spec = spheres[sphr].spec;
        diffuse = spheres[sphr].diffuse;
        reflCoeff = spheres[sphr].refl;
        refract_coeff=spheres[sphr].refract_coeff;
        shine= spheres[sphr].shine;
        normal = spNormal;
        newObj = ObjectID(sphr,SPHERE);
    }

    Ray rflct = calculateReflection(ri,normal,incidentVector);

    Ray refract = calculateRefraction(ri,normal,incidentVector);


    Color ret(0,0,0);

    Color ambientLight = sourcePower.dotProduct(colr);
    Color diffuseLight = sourcePower.dotProduct(colr);
    Color specLight = sourcePower.dotProduct(colr);

    ret = ambientLight.scale(amb);

    s=lights.size();
    for(int i = 0; i < s; i++)
    {
        Ray lightRay(lights[i],incidentVector);



        if(hasObstacle(lights[i], incidentVector)==true)    continue;

        Ray lightRflct = calculateReflection(lightRay,normal,incidentVector);

        Ray lightRefract = calculateRefraction(lightRay,normal,incidentVector);


        double t1=lightRflct.dir.dotProduct(normal);
        double t2=ri.dir.scale(-1).dotProduct(lightRflct.dir);

        double t1_refract=lightRefract.dir.dotProduct(normal.scale(-1));
        double t2_refract=ri.dir.dotProduct(lightRefract.dir);


        double normal_ab=sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z);
        double rflct_ab=sqrt(lightRflct.dir.x*lightRflct.dir.x+lightRflct.dir.y*lightRflct.dir.y+lightRflct.dir.z*lightRflct.dir.z);
        double ri_ab=sqrt(ri.dir.x*ri.dir.x+ri.dir.y*ri.dir.y+ri.dir.z*ri.dir.z);


        double normal_refract_ab=sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z);
        double refract_ab=sqrt(lightRefract.dir.x*lightRefract.dir.x+lightRefract.dir.y*lightRefract.dir.y+lightRefract.dir.z*lightRefract.dir.z);
        //double ri_ab=sqrt(ri.dir.x*ri.dir.x+ri.dir.y*ri.dir.y+ri.dir.z*ri.dir.z);

        double cosTheta = max(0.0,t1/(rflct_ab*normal_ab));
        double cosPhi = max(0.0,t2/(ri_ab*rflct_ab));

        double refract_cosTheta = max(0.0,t1_refract/(refract_ab*normal_refract_ab));
        double refract_cosPhi = max(0.0,t2_refract/(ri_ab*refract_ab));




        double mul=diffuse*cosTheta;
        Color df = diffuseLight.scale(mul);


        mul=spec*pow(cosPhi,shine);
        Color spc = specLight.scale(mul);

        Color tmp=df.plusColor(spc);
        ret = ret.plusColor(tmp);

        double mul_refract=diffuse*refract_cosTheta;
        df = diffuseLight.scale(mul_refract);

        mul_refract=spec*pow(refract_cosPhi,shine);
        spc = specLight.scale(mul_refract);

        tmp=df.plusColor(spc);
        ret = ret.plusColor(tmp);
    }

    Color lght = rayCast(rflct,level+1, newObj);

    Color rfrct= rayCast(refract,level+1,newObj);

    ret = ret.plusColor(lght.scale(reflCoeff));


    ret = ret.plusColor(rfrct.scale(refract_coeff));


    ret.r = min(ret.r,1.0);
    ret.g = min(ret.g,1.0);
    ret.b = min(ret.b,1.0);
    return ret;
}

void rayTracer()
{
    Ray src;
    src.p0 = cameraPos;

    for(int row = 0; row<screenHeight; row++)
    {
        double ht = row-screenHeight/2;
        ht/=screenHeight/2;
        for(int colm = 0; colm<screenWidth; colm++)
        {
            double wd = colm-screenWidth/2;
            wd/=screenWidth/2;

            Vector t2=(cameraRight.scale(wd));
            Vector t1=(cameraUp.scale(ht));
            Vector temp = t1.plusVector(t2) ;
            Vector dir = cameraLook.plusVector(temp);

            src.dir = dir;
            imageMap[colm][row] = rayCast(src,0, ObjectID(0,EYE));

        }
    }

}

void display()
{

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(cameraPos.x,cameraPos.y,cameraPos.z,
              cameraPos.x+100*cameraLook.x,cameraPos.y+100*cameraLook.y,cameraPos.z+100*cameraLook.z,
              cameraUp.x,cameraUp.y,cameraUp.z);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    {

        drawAxes(1000.0);

        int s=spheres.size();
        for(int i = 0; i < s; i++)
        {
            spheres[i].draw();
        }

        s=triangles.size();
        for(int i = 0; i < s; i++)
        {
            triangles[i].draw();
        }

        s=lights.size();
        for(int i = 0; i <s; i++)
        {
            drawCube(lights[i]);
        }

        checkBoard.draw();

    }
    glPopMatrix();
    glutSwapBuffers();
}

void animate()
{
    glutPostRedisplay();
}

void init()
{
    cameraPos=Vector(0,-150,10);
    cameraLook=Vector(0,1,0);
    cameraUp=Vector(0,0,1);
    cameraRight=Vector(1,0,0);

    glClearColor(0,0,0,0);


    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
    gluPerspective(85,	1,	1,	1000.0);

}

void readInput()
{
    int totalObj;
    string type;

    scanf("%lf %lf",&recursionLevel,&screenWidth);
    scanf("%d",&totalObj);

    screenHeight=screenWidth;
    checkBoard = CheckBoard(screenHeight,screenWidth);

    sourcePower.r=1.0;
    sourcePower.g=1.0;
    sourcePower.b=1.0;

    texture=false;

    imageMap= new Color*[(int)floor(screenWidth)+2];
    checkImage= new Color*[(int)floor(screenWidth)+2];

    for(int i=0; i<floor(screenWidth)+2; i++)
    {

        imageMap[i]=new Color[(int)floor(screenHeight)+2];
        checkImage[i]=new Color[(int)floor(screenHeight)+2];

    }


    for(int xx = 0;xx<(int)floor(screenWidth);xx++)
    {
        for(int yy = 0;yy<(int)floor(screenHeight);yy++)
        {
            int row=xx/checkBoard.boardWidth;
            int colm=yy/checkBoard.boardHeight;

            if((row+colm)%2)
            {
                double c=1.0;
                checkImage[xx][yy]=Color(c,c,c);
            }
            else
            {
                double c=0.0;
                checkImage[xx][yy]=Color(c,c,c);
            }
        }
    }


    for(int i = 0; i < totalObj; i++)
    {
        cin>>type;

        Color c;
        double am,df,spc, refl,shn;

        if(type=="sphere")
        {
            Sphere sp;
            sp.read();
            spheres.push_back(sp);
        }
        else if(type == "general")
        {
            double A,B,C,D,E,F,G,H,I,J;
            double x,y,z;
            double l,w,h;

            scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&A,&B,&C,&D,&E,&F,&G,&H,&I,&J);
            scanf("%lf %lf %lf",&x,&y,&z);
            scanf("%lf %lf %lf",&l,&w,&h);
            scanf("%lf %lf %lf",&c.r,&c.g,&c.b);
            scanf("%lf %lf %lf %lf %lf",&am,&df,&spc,&refl,&shn);


            General gn(A,B,C,D,E,F,G,H,I,J);
            gn.cube_reference=Vector(x,y,z);
            gn.l=l;
            gn.w=w;
            gn.h=h;
            gn.c=c;
            gn.amb = am;
            gn.diffuse = df;
            gn.spec = spc;
            gn.refl = refl;
            gn.shine = shn;
            gn.refract_coeff=min(refl*eff,1.0);
            generals.push_back(gn);
        }
        else if(type=="triangle")
        {
            Vector a1,b1,c1;


            scanf ( "%lf%lf%lf",&a1.x,&a1.y,&a1.z);
            scanf ( "%lf%lf%lf",&b1.x,&b1.y,&b1.z);
            scanf ( "%lf%lf%lf",&c1.x,&c1.y,&c1.z);

            scanf("%lf %lf %lf",&c.r,&c.g,&c.b);

            scanf("%lf %lf %lf %lf %lf",&am,&df,&spc,&refl,&shn);



            Triangle t=Triangle(a1,b1,c1);

            t.c=c;
            t.amb = am;
            t.diffuse = df;
            t.spec = spc;
            t.refl = refl;
            t.shine = shn;
            t.refract_coeff=min(refl*eff,1.0);
            triangles.push_back(t);
        }

    }
    int lightCnt;
    scanf("%d",&lightCnt);
    for(int i = 0; i < lightCnt; i++)
    {

        double xx,yy,zz;
        Vector lg;
        cin>>xx>>yy>>zz;
        lg=Vector(xx,yy,zz);
        lights.push_back(lg);

    }

}

void printImage()
{

    bitmap_image image((int)floor(screenWidth),(int)floor(screenHeight));

    for(int i=0; i<(int)floor(screenWidth); i++)
    {
        for(int j=0; j<(int)floor(screenHeight); j++)
        {
            image.set_pixel(i,j,imageMap[i][(int)floor(screenHeight)-1-j].r*255,imageMap[i][(int)floor(screenHeight)-1-j].g*255,imageMap[i][(int)floor(screenHeight)-1-j].b*255);
        }
    }

    image.save_image("out.bmp");

}

int main(int argc, char **argv)
{

    freopen("scene.txt","r",stdin);
    readInput();

    glutInit(&argc,argv);
    glutInitWindowSize(screenHeight,screenWidth);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

    glutCreateWindow("Ray-Tracer");

    init();

    glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(display);
    glutIdleFunc(animate);

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    glutMainLoop();

    return 0;
}

void drawAxes(double l)
{
    glBegin(GL_LINES);
    {
        glColor3f(0, 0, 1.0);
        glVertex3f( l,0,0);
        glVertex3f(-l,0,0);

        glColor3f(0, 1.0, 0);
        glVertex3f(0,-l,0);
        glVertex3f(0, l,0);

        glColor3f(1.0, 0, 0);
        glVertex3f(0,0,-l);
        glVertex3f(0,0, l);

    }
    glEnd();


}

void drawSquare(double x,double y,double z,double a,double c)
{
    glPushMatrix();
    {
        glTranslatef(x,y,z);
        glRotatef(90,0,1,0);
        a/=2;
        glColor3f(c,c,c);
        glBegin(GL_QUADS);
        {
            glVertex3f( 0,a, a);
            glVertex3f( 0,a,-a);
            glVertex3f( 0,-a,-a);
            glVertex3f( 0,-a, a);
        }
        glEnd();
    }
    glPopMatrix();

}

void drawSquare1(double x,double y,double z,double a,double c)
{
    glPushMatrix();
    {
        glTranslatef(x,y,z);
        glRotatef(90,0,1,0);
        a/=2;
        glColor3f(0,c,0);
        glBegin(GL_QUADS);
        {
            glVertex3f( 0,a, a);
            glVertex3f( 0,a,-a);
            glVertex3f( 0,-a,-a);
            glVertex3f( 0,-a, a);
        }
        glEnd();
    }
    glPopMatrix();

}

void drawCube(Vector centre)
{
    glPushMatrix();
    {
        glTranslatef(centre.x,centre.y,centre.z);
        drawSquare1(0,0,-0.5,1.0,1.0);
        drawSquare1(0,0,0.5,1.0,1.0);

        glPushMatrix();
        {
            glTranslatef(0,-0.5,0);
            glRotatef(90,1,0,0);
            drawSquare1(0,0,0,1.0,1.0);
        }
        glPopMatrix();

        glPushMatrix();
        {
            glTranslatef(0,0.5,0);
            glRotatef(90,1,0,0);
            drawSquare1(0,0,0,1.0,1.0);
        }
        glPopMatrix();

        glPushMatrix();
        {
            glTranslatef(0.5,0,0);
            glRotatef(90,0,1,0);
            drawSquare1(0,0,0,1.0,1.0);
        }
        glPopMatrix();

        glPushMatrix();
        {
            glTranslatef(-0.5,0,0);
            glRotatef(90,0,1,0);
            drawSquare1(0,0,0,1.0,1.0);
        }
        glPopMatrix();
    }
    glPopMatrix();

}

void drawSphere(double radius,int slices,int stacks,double x,double y,double z)
{
    glPushMatrix();
    {

        glTranslatef(x,y,z);

        stacks*=4;
        slices*=4;
        Vector points[stacks+2][slices+2];
        int i,j;
        double h,r;
        for(i=0; i<=stacks; i++)
        {
            double ang=((double)(stacks-i)/(double)stacks)*(pi*2);
            r=radius*cos(-ang);
            h=radius*sin(-ang);
            for(j=0; j<=slices; j++)
            {
                points[i][j].z=h;
                double ang2=((double)(slices-j)/(double)slices)*pi*2;
                points[i][j].y=r*sin(-ang2);
                points[i][j].x=r*cos(-ang2);
            }
        }

        for(i=0; i<stacks; i++)
        {
            for(j=0; j<slices; j++)
            {
                glBegin(GL_TRIANGLES);
                {

                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);

                }
                glEnd();

                glBegin(GL_TRIANGLES);
                {
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                }
                glEnd();
            }
        }
    }
    glPopMatrix();
}

void keyboardListener(unsigned char key, int x,int y)
{
    Vector l=cameraLook,r=cameraRight,u=cameraUp;
    double ang=.05;

    switch(key)
    {
    case '1': // object right

        cameraLook.x=l.x*cos(ang)-r.x*sin(ang);
        cameraLook.y=l.y*cos(ang)-r.y*sin(ang);
        cameraLook.z=l.z*cos(ang)-r.z*sin(ang);

        cameraRight=cameraLook.crossVector(cameraUp);


        break;
    case '2':  //object left
        ang*=-1;
        cameraLook.x=l.x*cos(ang)-r.x*sin(ang);
        cameraLook.y=l.y*cos(ang)-r.y*sin(ang);
        cameraLook.z=l.z*cos(ang)-r.z*sin(ang);

        cameraRight=cameraLook.crossVector(cameraUp);
        break;
    case '3':  //object up
        cameraUp.x=u.x*cos(ang)-l.x*sin(ang);
        cameraUp.y=u.y*cos(ang)-l.y*sin(ang);
        cameraUp.z=u.z*cos(ang)-l.z*sin(ang);

        cameraLook=cameraUp.crossVector(cameraRight);
        break;
    case '4':  //object down
        ang*=-1;
        cameraUp.x=u.x*cos(ang)-l.x*sin(ang);
        cameraUp.y=u.y*cos(ang)-l.y*sin(ang);
        cameraUp.z=u.z*cos(ang)-l.z*sin(ang);

        cameraLook=cameraUp.crossVector(cameraRight);
        break;

    case '5':  // clockwise
        cameraRight.x=r.x*cos(ang)-u.x*sin(ang);
        cameraRight.y=r.y*cos(ang)-u.y*sin(ang);
        cameraRight.z=r.z*cos(ang)-u.z*sin(ang);

        cameraUp=cameraRight.crossVector(cameraLook);
        break;
    case '6': //anticlockwise
        ang*=-1;
        cameraRight.x=r.x*cos(ang)-u.x*sin(ang);
        cameraRight.y=r.y*cos(ang)-u.y*sin(ang);
        cameraRight.z=r.z*cos(ang)-u.z*sin(ang);

        cameraUp=cameraRight.crossVector(cameraLook);
        break;
    case '0':
        rayTracer();
        printImage();
        break;
    case '7':
        changeBoard();
        rayTracer();
        printImage();
        break;
    }




}


void specialKeyListener(int key, int x,int y)
{
    double mov=3;
    Vector ret;
    switch(key)
    {

    case GLUT_KEY_PAGE_DOWN:		//down arrow key

        ret=cameraPos.minusVector(cameraUp.scale(mov));
        cameraPos=Vector(ret.x,ret.y,ret.z);

        //cout<<cameraPos.x<<" "<<cameraPos.y<<" "<<cameraPos.z<<endl;
        break;
    case GLUT_KEY_PAGE_UP:		// up arrow key
        ret=cameraPos.plusVector(cameraUp.scale(mov));
        cameraPos=Vector(ret.x,ret.y,ret.z);

        //cout<<cameraPos.x<<" "<<cameraPos.y<<" "<<cameraPos.z<<endl;

        break;

    case GLUT_KEY_RIGHT:

        ret=cameraPos.plusVector(cameraRight.scale(mov));
        cameraPos=Vector(ret.x,ret.y,ret.z);

        //cout<<cameraPos.x<<" "<<cameraPos.y<<" "<<cameraPos.z<<endl;

        break;
    case GLUT_KEY_LEFT:
        ret=cameraPos.minusVector(cameraRight.scale(mov));
        cameraPos=Vector(ret.x,ret.y,ret.z);

        //cout<<cameraPos.x<<" "<<cameraPos.y<<" "<<cameraPos.z<<endl;

        break;

    case GLUT_KEY_UP:
        ret=cameraPos.plusVector(cameraLook.scale(mov));
        cameraPos=Vector(ret.x,ret.y,ret.z);

        //cout<<cameraPos.x<<" "<<cameraPos.y<<" "<<cameraPos.z<<endl;

        break;
    case GLUT_KEY_DOWN:
        ret=cameraPos.minusVector(cameraLook.scale(mov));
        cameraPos=Vector(ret.x,ret.y,ret.z);

        //cout<<cameraPos.x<<" "<<cameraPos.y<<" "<<cameraPos.z<<endl;

        break;
    }

}


