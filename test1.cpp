#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

static inline double sqr(double x) { return x * x; }

//using namespace std;

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
    double& operator[](int i) { return coord[i]; }
    double operator[](int i) const { return coord[i]; }

    Vector& operator+=(const Vector& v) {
        coord[0] += v[0];
        coord[1] += v[1];
        coord[2] += v[2];
        return *this;
    }

    double norm2() const {
        return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
    }

    void normalize()
    {
        double norme = sqrt(norm2());
        coord[0] /= norme;
        coord[1] /= norme;
        coord[2] /= norme;
    }

    double coord[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
 
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class Ray
{
    public:

        Ray(const Vector& O,const Vector& u) : O(O),u(u)
        {};

        Vector O;
        Vector u;


};

class Sphere 
{
    public:
        
        Sphere(const Vector& C,double R, Vector albedo) : C(C),R(R),albedo(albedo)
        {};

        bool intersect(const Ray& r, Vector& P, Vector& N)
        {
            // solve at^2 + bt +c =0
            double a = 1;
            double b = 2 * dot(r.u,r.O-C);
            double c = (r.O-C).norm2()-R*R;

            double delta = b*b-4*a*c;

            if (delta <= 0) return false;

            double t;
            
            double sqrtDelta = sqrt(delta);
            double t0 = (-b -sqrtDelta)/2;
            double t1 = (-b +sqrtDelta)/2;

            if (t1 < 0) return false;
            if (t0 > 0) 
            {
                t = t0;
            }
            else{
                t = t1;
            }
            P = r.O + t*r.u;
            N = (P-C);
            N.normalize();
            return true;


        };

        Vector C;
        double R;
        Vector albedo;
};

class Scene
{
    public :
        Scene() {}
        void addSphere(const Sphere& s) {objects.push_back(s);};

        bool intersect(const Ray& r, Vector& P, Vector& N, int& id)
        {
            bool has_inter = false;

            for (int i = 0; i < objects.size(); i++)
            {
                Vector local_N, local_P;
                if (objects[i].intersect(r, local_P, local_N))
                {
                    
                    if (has_inter == false || ((P-r.O).norm2() >= (local_P-r.O).norm2()))
                    {
                        N  = local_N;
                        P  = local_P;
                        id = i;
                    }
                    
                    has_inter = true;

                }
            }
            
            return has_inter;
        }


        std::vector<Sphere> objects;
        
};

int main() {
    int W = 512;
    int H = 512;
    
    Scene scene;

    Sphere sph1(Vector (0,0,0),10,Vector (1.,1.,1.));
    Sphere  plafond(Vector (0., 1000.,0.),1000-60,Vector (0.3,0.3,0.3));
    Sphere      sol(Vector (0.,-1000.,0.),1000-10,Vector (0.3,0.3,0.3));
    Sphere    mur_G(Vector ( 1000.,0.,0.),1000-60,Vector (1.,0.,0.));
    Sphere    mur_D(Vector (-1000.,0.,0.),1000-60,Vector (0.,0.,1.));
    Sphere derriere(Vector (0.,0., 1000.),1000-60,Vector (0.5,0.,0.5));
    Sphere   devant(Vector (0.,0.,-1000.),1000-60,Vector (0.,1.,0.));

    scene.addSphere(sph1);
    scene.addSphere(plafond);
    scene.addSphere(sol);
    scene.addSphere(mur_G);
    scene.addSphere(mur_D);
    scene.addSphere(derriere);
    scene.addSphere(devant);

    
    const Vector source(-10,20,40);
    
    const Vector camera(0, 0, 55);
    const double fov = 60 *M_PI /180;
    const double I = 1E10;


 
    std::vector<unsigned char> image(W*H * 3, 0);
    

#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            
            Vector P, N;
            int id;

            Vector u(j-W/2+0.5, H/2-i-0.5, -W/(2.*tan(fov/2)));
            u.normalize();

            Ray R(camera,u);

            if (scene.intersect(R,P,N,id))
            {
                Vector veclum = source-P;
                float normVeclum = veclum.norm2();
                veclum.normalize();
                Vector color = ( (std::max(0.,dot(N,veclum))*I / ((M_PI*4*M_PI*normVeclum)) ) *scene.objects[id].albedo);
                image[(i*W + j) * 3 + 0] = std::min(255.,pow(color[0],0.45));    // RED
                image[(i*W + j) * 3 + 1] = std::min(255.,pow(color[1],0.45));   // GREEN
                image[(i*W + j) * 3 + 2] = std::min(255.,pow(color[2],0.45));  // BLUE
            }
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}