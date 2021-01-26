#ifndef _HEADER_H_
#define _HEADER_H_

#include <cmath>
#include <cstdio>
#include <vector>
#include <time.h>
#include <sys/stat.h>
#include <iostream>

// 3D ベクトルの class
struct Vec3D{
  double x;
  double y;
  double z;
  // デフォルトコンストラクタ
  Vec3D() = default;

  //Vec3D v{10,20}; などと初期化できるように
  Vec3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z){}

  // const:constオブジェクトに大しても実行できるように
  double length() const {
    return sqrt(x*x+y*y+z*z);
  }

  Vec3D operator +() const{
    return *this;
  }
  Vec3D operator -() const{
    return {-x,-y,-z};
  }
  Vec3D operator +(const Vec3D &v) const{
    return {x+v.x, y+v.y, z+v.z};
  }
  Vec3D operator -(const Vec3D &v) const{
    return {x-v.x, y-v.y, z-v.z};
  }
  Vec3D operator *(double s) const{
    return {x*s, y*s, z*s};
  }
  Vec3D operator /(double s) const{
    return {x/s, y/s, z/s};
  }
};

inline Vec3D operator *(double s, const Vec3D &v){
  return {s*v.x, s*v.y, s*v.z};
}

//外積 Vec3D%Vec3D
inline Vec3D operator%(const Vec3D& u,const Vec3D& v){
	Vec3D w;
	w.x=u.y*v.z-u.z*v.y;
	w.y=u.z*v.x-u.x*v.z;
	w.z=u.x*v.y-u.y*v.x;
	return w;
}

class Particle_list{
protected:
  class Particle;
  char FILE_PATH[100];
  void init_file();
  void init_particles();
public:
  std::vector<Particle> P_list;
  Particle_list(int n){
    P_list = std::vector<Particle>(n);
    init_file();
    init_particles();
  }
  void output(double t);
};

class Particle_list::Particle{
public:
  Vec3D r;
  Vec3D v;
  double Energy() const{
    return 0.5*v.length()*v.length()-1.0/r.length();    
  }
  Particle() = default;
  void init(double r0, double phi);
};

// z 方向単位ベクトル
const Vec3D e_z = {0,0,1};

Vec3D force(Vec3D const &r, Vec3D const &v);
void update(Particle_list &P_list, double &t, double dt);
void euler(Vec3D &r, Vec3D &v, double dt);
void runge_kutta(Vec3D &r, Vec3D &v, double dt);
void calc_time(int t);
double cot(double x);

#endif