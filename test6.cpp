#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


#include <string>
#include <stdio.h>
#include <algorithm>

#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include <random>
#include <omp.h>

using namespace std;

static inline double sqr(double x) { return x * x; }

std::default_random_engine engine[8];
std::uniform_real_distribution<double> uniform(0.0,1.0);

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
Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const Vector& a, double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator*(double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
 
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]);
}

Vector random_cos(const Vector& N)
{
    int thread_id = omp_get_thread_num();
    double r1 = uniform(engine[thread_id]);
    double r2 = uniform(engine[thread_id]);

    double r= sqrt(1-r2);
    double x = cos(2.*M_PI * r1) * r;
    double y = sin(2.*M_PI * r1) * r;
    double z = sqrt(r2);

    Vector T1;
    if (abs(N[0]) < abs(N[1]) && abs(N[0]) < abs(N[2]))
    {
        T1 = Vector(0, -N[2], N[1]);
    }
    else if (abs(N[1]) < abs(N[0]) && abs(N[1]) < abs(N[2]))
    {
        T1 = Vector(-N[2], 0 , N[0]);
    }
    else  
    {
        T1 = Vector(-N[1], N[0], 0);
    }

    T1.normalize();
    Vector T2 = cross(N,T1);

    return z*N + x*T1 + y*T2;
}

class Ray
{
    public:

        Ray(const Vector& O,const Vector& u) : O(O),u(u)
        {};
      

        Vector O;
        Vector u;


};

class Object
{
    public:

        Object()  {};
        Object(Vector albedo, bool mirror, bool transparent, double refraction, bool creuse) : albedo(albedo), mirror(mirror), transparent(transparent), refraction(refraction), creuse(creuse)  {};

        virtual bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const = 0;

        Vector albedo;
        bool mirror;
        bool transparent;
        double refraction;
        bool creuse;
};

class Sphere : public Object
{
    public:
        
        Sphere(const Vector& C,double R, Vector albedo, bool mirror, bool transparent, double refraction, bool creuse) : C(C),R(R), ::Object(albedo, mirror, transparent, refraction, creuse)
        {};
        Sphere(const Vector& C,double R, Vector albedo, bool mirror, bool transparent, double refraction) :              C(C),R(R), ::Object(albedo, mirror, transparent, refraction,  false)
        {};
        Sphere(const Vector& C,double R, Vector albedo) :                                                                C(C),R(R), ::Object(albedo,  false,       false,        1.0,  false)
        {};

        bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const
        {
            // solve at^2 + bt +c =0
            double a = 1;
            double b = 2 * dot(r.u,r.O-C);
            double c = (r.O-C).norm2()-R*R;

            double delta = b*b-4*a*c;

            if (delta <= 0) return false;
            
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
            if (creuse) 
            {
                N = -1*N;
            }

            return true;


        };

        Vector C;
        double R;
        
};

class Scene
{
    public :
        Scene() {}
        void addObject(const Object* s) {objects.push_back(s);};

        void addSource(const Vector& source, float i)
        {
            S = source;
            I = i;
        };

        bool intersect(const Ray& r, Vector& P, Vector& N, int& id, double& t) const
        {
            bool has_inter = false;

            for (int i = 0; i < objects.size(); i++)
            {
                Vector local_N, local_P;
                double local_t; 
                if (objects[i]->intersect(r, local_P, local_N, local_t))
                {
                    
                    if (has_inter == false || ((P-r.O).norm2() >= (local_P-r.O).norm2()))
                    {
                        N  = local_N;
                        P  = local_P;
                        t  = local_t;
                        id = i;
                    }
                    
                    has_inter = true;

                }
            }
            
            return has_inter;
            }

        Vector getColor(const Ray& ray)
        {
            return getColor(ray,20);
        }

        Vector getColor(const Ray& R,const int n_rebond)
        {
            Vector noir(0,0,0);

            if (n_rebond <= 0) return noir;

            Vector N, P;
            int id;
            double t;
            
            Vector lumC = dynamic_cast<const Sphere*> (this->objects[0])->C;
            double lumR = dynamic_cast<const Sphere*> (this->objects[0])->R;
            
            if (intersect(R, P, N, id, t))
            {
                if (id == 0)
                {
                    double intensite = I/(4*sqr(M_PI*lumR));
                    return Vector(intensite,intensite,intensite);
                }
                else if (objects[id]->mirror)
                {
                    Vector reflected_vector = R.u - 2*dot(N,R.u)*N;
                    Ray reflected_ray(P + 0.001*N,reflected_vector);

                    return getColor(reflected_ray, n_rebond-1);
                }
                else if (objects[id]->transparent)
                {
                    double n1 = refraction_indice;
                    double n2 = refraction_indice;
                    Vector N_correct;

                    if (dot(N, R.u) < 0 )
                    {
                        n2 = objects[id]->refraction;
                        N_correct = N;
                    }
                    else
                    {
                        n1 = objects[id]->refraction;
                        N_correct = -1*N;
                        
                    }
                    
                    Vector WtT = (n1/n2)*(R.u - dot(R.u,N_correct)*N_correct);
                    Vector WtN(0,0,0);
                    if (1-sqr(n1/n2)*(1-sqr(dot(R.u,N_correct))) > 0)
                    {
                        WtN = -1*sqrt(1-sqr(n1/n2)*(1-sqr(dot(R.u,N_correct))))*N_correct;
                    }
                     

                    Ray refracted_ray(P-0.001*N_correct,WtT+WtN);
                    return getColor(refracted_ray, n_rebond-1);

                }
                else
                {
                    // eclairage direct

                    /*Vector veclum = S-P;
                    float normVeclum = veclum.norm2();
                    veclum.normalize();

                    Vector color  = noir;
                    
                    int idprime;
                    Vector Pprime, Nprime;
                    double tprime;

                    bool has_intersect = intersect(Ray(P+0.001*N, veclum), Pprime, Nprime, idprime, tprime);
                    bool shadow = has_intersect && (tprime*tprime < normVeclum);

                    if (!shadow) 
                    {
                        color  = ( (std::max(0.,dot(N,veclum))*I / ((M_PI*4*M_PI*normVeclum)) ) *objects[id].albedo);
                    }
                    */

                    // eclairage direct, ombre douce
                   
                    Vector dirLum = lumC - P;
                    dirLum.normalize();
                    Vector dir_xprime = random_cos(-dirLum);
                    Vector xprime = dir_xprime*lumR + lumC;

                    double proba_xprime = dot(-dirLum, dir_xprime)/(M_PI * sqr(lumR));

                    Vector xxprime = xprime - P;
                    double d2lum = xxprime.norm2();
                    xxprime.normalize();
                    


                    int idprime;
                    Vector Pprime, Nprime;
                    double tprime;

                    bool interlum = intersect(Ray(P+0.001*N, xxprime), Pprime, Nprime, idprime, tprime);
                    bool shadow = interlum && ( sqr(tprime) < d2lum)*0.95;
                    
                    double weighted_intensity = I/(4*sqr(M_PI*lumR));

                    Vector color  = noir;

                    

                    if (!shadow) 
                    {
                        //color  = (max(0.,dot(N,xxprime))*(max(0.,dot(xprime,-xxprime)))*weighted_intensity / (d2lum*proba_xprime*M_PI) )*objects[id].albedo;
                        color  = (1*1*weighted_intensity / (d2lum*proba_xprime*M_PI) )*objects[id]->albedo;
                    }

                    // eclairage indirect
                    Vector indirectColor= noir;
                    
                    Vector omega_i = random_cos(N);
                    indirectColor = objects[id]->albedo*getColor(Ray(P+0.001*N, omega_i),n_rebond-1); 
                    
                    color += indirectColor;
                    
                    return color;
                }
    
            }
            else
            {
                return noir;
            }
        }


    std::vector<const Object*> objects;
    Vector S;
    double I = 1E10;
    double refraction_indice = 1;
        
};

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};
 
class BoundingBox 
{
    public :
        bool intersect(const Ray& R) const
        {
            double tx1 = (bmin[0]-R.O[0]) / R.u[0];
            double tx2 = (bmax[0]-R.O[0]) / R.u[0];
            double txmin = min(tx1,tx2);
            double txmax = max(tx1,tx2);

            double ty1 = (bmin[1]-R.O[1]) / R.u[1];
            double ty2 = (bmax[1]-R.O[1]) / R.u[1];
            double tymin = min(ty1,ty2);
            double tymax = max(ty1,ty2);

            double tz1 = (bmin[2]-R.O[2]) / R.u[2];
            double tz2 = (bmax[2]-R.O[2]) / R.u[2];
            double tzmin = min(tz1,tz2);
            double tzmax = max(tz1,tz2);

            return txmin < tymax && txmin < tzmax && tymin < txmax && tymin < tzmax && tzmin < tymax && tzmin < txmax;
        }
    
        Vector bmin;
        Vector bmax;
};

class BVH
{
    public:

        BoundingBox bbox;
        int debut;
        int fin;
        BVH* filsg;
        BVH* filsd;

};

class TriangleMesh : public Object
{
public:
  ~TriangleMesh() {}

    TriangleMesh(Vector albedo, bool mirror, bool transparent, double refraction, bool creuse) : ::Object(albedo, mirror, transparent, refraction, creuse)  
    {};
    TriangleMesh(Vector albedo, bool mirror, bool transparent, double refraction) :              ::Object(albedo, mirror, transparent, refraction,  false)
    {};
    TriangleMesh(Vector albedo) :                                                                ::Object(albedo,  false,       false,        1.0,  false)
    {};

    void transform(double scale, Vector translation)
    {
        int len = vertices.size();

        for (int i = 0; i < len; i++)
        {
            vertices[i] = (vertices[i] * scale) + translation;
        }
    }
    
    void readOBJ(const char* obj) 
    { 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);
    }
    
    void compute_bbox(int debut, int fin)
    {
        bbox.bmax = Vector( -1E9, -1E9, -1E9);
        bbox.bmin = Vector(  1E9,  1E9, 1E9);
        int len = vertices.size();
        for (int i = debut; i < fin; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (bbox.bmax[j] < vertices[i][j]) {bbox.bmax[j] = vertices[i][j];} 
                if (bbox.bmin[j] > vertices[i][j]) {bbox.bmin[j] = vertices[i][j];} 
            }
        }
    }

    void compute_bbox()
    {
        bbox.bmax = Vector( -1E9, -1E9, -1E9);
        bbox.bmin = Vector(  1E9,  1E9, 1E9);
        int len = vertices.size();
        for (int i = 0; i < len; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (bbox.bmax[j] < vertices[i][j]) {bbox.bmax[j] = vertices[i][j];} 
                if (bbox.bmin[j] > vertices[i][j]) {bbox.bmin[j] = vertices[i][j];} 
            }
        }
    }

 /*   intersection:
        std::list<Noeud*> l;
        l.push(racine)
        tant que l est pas vide
            noeud n = l.pop(
                if (n-> filsg)
                {

                }
            )*/
/*
    void calculBVH(BVH* node, int debut, int fin)
    {
        n->indices = 
        n->boite =compute_bbox

        axis = + grand axe

        pivot = debut
        pour chaque triangle entre début et fin
            si milieu triangle < milieu boite
                swap(triangle, pivot)
                pivot++;
            si pivot = debut et fin 
                noeud*g,d,n -> filsg = g, filsd = d
                calculBVH(g,debut,pivot)
                calculBVH(d,pivot,fin)

    }*/

    bool intersect(const Ray& R, Vector& P, Vector& N, double& t) const
    {
        if (!bbox.intersect(R)) return false;

        bool has_inter = false;
        int len = indices.size();
        for (int i = 0; i < len; i++)
        {
            Vector A = vertices[indices[i].vtxi];
            Vector B = vertices[indices[i].vtxj];
            Vector C = vertices[indices[i].vtxk];

            Vector e1 = B-A;
            Vector e2 = C-A;

            Vector local_N = cross(e1,e2);

            //if (dot(local_N,R.u) > 0) local_N = -1 * local_N; 

            double denom = dot(R.u, local_N);

            double local_t = dot(A-R.O, local_N)/denom;
            if (local_t < 0) continue;

            double beta = -dot(R.u, cross(A-R.O, e2)) / denom;
            if (beta < 0 || beta > 1) continue;

            double gamma = dot(R.u, cross(A-R.O, e1)) / denom;
            if (gamma < 0 || gamma > 1) continue;

            double alpha = 1 - beta - gamma;
            if (alpha < 0) continue;

            if (local_t < t)
            {
                t = local_t;
                has_inter = true;
                N = local_N;
            }
        }       

        N.normalize();
        P = R.O + t*R.u;
        return has_inter;
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;

    BoundingBox bbox = BoundingBox();
    
};

    

int main() {
    int W = 512;
    int H = 512;

    Scene scene;

    Sphere  source(Vector ( -10, 20 , 40),     10,Vector (0.7,0.7,0.7));

    Sphere     sph2(Vector (0  ,  0 ,  0),    9.5,Vector (0.7,0.7,0.7) ,false ,   true,1.3,true);
    Sphere     sph3(Vector (-20 ,  0 ,  0),    10,Vector (0.7,0.7,0.7) ,true  ,  false,  1);
    Sphere     sph4(Vector ( 20 ,  0 ,  0),    10,Vector (0.7,0.7,0.7) ,true  ,  false,  1);

    Sphere     sph1(Vector (0  ,  0 ,  0),     10,Vector (0.7,0.7,0.7), false,  true, 1.3);

    Sphere  plafond(Vector (0., 1000.,0.),1000-60,Vector (0.3,0.3,0.3));
    Sphere      sol(Vector (0.,-1000.,0.),1000-10,Vector (0.3,0.3,0.3));
    Sphere    mur_G(Vector ( 1000.,0.,0.),1000-60,Vector (1. , 0.,0. ));
    Sphere    mur_D(Vector (-1000.,0.,0.),1000-60,Vector (0. ,0. ,1. ));
    Sphere derriere(Vector (0.,0., 1000.),1000-60,Vector (0.5,0. ,0.5));
    Sphere   devant(Vector (0.,0.,-1000.),1000-60,Vector (0. ,1.,0.  ));

    scene.addObject(&source);

    //scene.addObject(&sph1);

    //scene.addObject(&sph2);
    //scene.addObject(&sph3);
    //scene.addObject(&sph4);

    scene.addObject(&plafond);
    scene.addObject(&sol);
    scene.addObject(&mur_G);
    scene.addObject(&mur_D);
    scene.addObject(&derriere);
    scene.addObject(&devant);

    TriangleMesh mesh(Vector(1.,1.,1.));
    mesh.readOBJ("cat.obj");
    mesh.transform(0.6,Vector(0.,-10.,0));
    mesh.compute_bbox();
    scene.addObject(&mesh);
 
    const Vector camera(0, 0, 55);
    const double fov = 60 *M_PI /180;

    const double dist_focus = 45;
 
    std::vector<unsigned char> image(W*H * 3, 0);

    int Nrays = 120;

#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        int tid = omp_get_thread_num(); 
        for (int j = 0; j < W; j++) {
            if (j==0) std::cout << i << std::endl;
            
            
            Vector P, N;
            int id;
            double t;

            Vector u;
            
            Vector color = Vector(0.,0.,0.);
            for (int k = 0; k < Nrays; k++)
            {   
                double r1 = uniform(engine[tid]);
                double r2 = uniform(engine[tid]);
                double r  = sqrt(-2.*log(r1));
                double gx = r*cos(2.*M_PI*r2)*1.;
                double gy = r*sin(2.*M_PI*r2)*1.;

                u = Vector(j-W/2.+0.5+gx, H/2.-i-0.5+gy, -W/(2.*tan(fov/2.)));
                u.normalize();

                Ray R(camera,u);
                color += scene.getColor(R);
            }
            color = color*(1.0/Nrays);
            

            image[(i*W + j) * 3 + 0] = std::min(255.,pow(color[0],0.45));    // RED
            image[(i*W + j) * 3 + 1] = std::min(255.,pow(color[1],0.45));   // GREEN
            image[(i*W + j) * 3 + 2] = std::min(255.,pow(color[2],0.45));  // BLUE
        }
    }
                
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}
