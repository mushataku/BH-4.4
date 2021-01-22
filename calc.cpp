#include "header.h"

/******************config********************/
// const int INITIAL = 0;
// 0:重力のみ 1:Lorentzt力のみ(半径依存無し) 2:Lorentzt力のみ(半径依存有) 3:重力+Lorentz
const int FORCE = 3;
// 0:Euler 1:Runge-Kutta 2:Leap Flog
const int METHOD = 1;
/******************config********************/

/******************計算条件********************/
const double beta = 0.0;
// const double beta = M_PI/10;
// スピンパラメータ
const double a_spin = 1.0;
const int N_particle = 50;
const double r0_ring = 10.0;
const int nt = 1000001;
const int dn = 10000;
double tmin = 0.0;
double tmax = 2.0*M_PI*std::pow(r0_ring, 1.5)*5.2;
// double tmax = M_PI/(2.0*a_spin)*10;
// double tmax = M_PI/(2.0*a_spin)*r0_ring*r0_ring*r0_ring*3;
double dt = (tmax-tmin)/double(nt-1);
/******************計算条件********************/

int main(){
  // 時間計測
  clock_t start_t, end_t;
  start_t = time(NULL);

  double t = tmin;
  Particle_list P_list(N_particle);
  output_data output;

  // 初期配置
  init_particle(P_list);
  // 出力
  output.write(P_list, t);

  printf("****************CALICULATION START****************\n");

  for(int i = 0; i < nt; i++) {
    update(P_list, t, dt);
    if(i%dn == 0) output.write(P_list, t);
  }

  printf("******************NORMAL END******************\n");
  end_t = time(NULL);
  // printf("This calculatioin took %ld second \n\n", end_t - start_t);
  output.calc_time((int)(end_t-start_t));
  return 0;
}

/******************関数********************/

void init_particle(Particle_list &P_list){
  // if(INITIAL == 0){
    double phi, r0, v0, x, y, z;
    for(int i = 0; i < N_particle; i++) {
      r0 = r0_ring;
      v0 = 1.0/std::sqrt(r0_ring);
      // v0 = 4.0*a_spin*r0;
      // v0 = 4.0*a_spin/(r0*r0);
      phi = 2.0*M_PI*i/N_particle;
      P_list[i].r = {r0*std::cos(phi), r0*std::cos(beta)*std::sin(phi), 
                          r0*std::sin(beta)*std::sin(phi)};
      P_list[i].v = {-v0*std::sin(phi), v0*std::cos(beta)*std::cos(phi), 
                          v0*std::sin(beta)*std::cos(phi)};
    // }
  }
}

// 重力だけ
Vec3D force(Vec3D const &r, Vec3D const &v){
  double len = r.length();
  Vec3D ret;
  if(FORCE == 0) ret = -r/(len*len*len);
  if(FORCE == 1){
    ret = -4.0*a_spin*v%e_z;
  }
  if(FORCE == 2){
    ret = -4.0*a_spin*v%e_z/(len*len*len);
  }
  if(FORCE == 3){
    ret = -(r + 4.0*a_spin*v%e_z)/(len*len*len);
  }
  return ret;
}

void update(Particle_list &P_list, double &t, double dt){
  for(int i = 0; i < N_particle; i++) {
    if(METHOD == 0) euler(P_list[i].r, P_list[i].v, dt);
    if(METHOD == 1) runge_kutta(P_list[i].r, P_list[i].v, dt);
  }
  t += dt;
}

void euler(Vec3D &r, Vec3D &v, double dt){
  Vec3D rk1,vk1;
  rk1 = v;
  vk1 = force(r, v);

  r = r + rk1*dt;
  v = v + vk1*dt;
}

void runge_kutta(Vec3D &r, Vec3D &v, double dt){
  Vec3D rk1,vk1,rk2,vk2, rk3,vk3,rk4,vk4;

  rk1 = v;
  vk1 = force(r, v);

  rk2 = v + 0.5*dt*vk1;
  vk2 = force(r + 0.5*dt*rk1, v + 0.5*dt*vk1);

  rk3 = v + 0.5*dt*vk2;
  vk3 = force(r + 0.5*dt*rk2, v + 0.5*dt*vk2);

  rk4 = v + dt*vk3;
  vk4 = force(r + dt*rk3, v + dt*vk3);

  r = r + dt*(rk1 + 2.0*rk2 + 2.0*rk3 + rk4)/6.0;
  v = v + dt*(vk1 + 2.0*vk2 + 2.0*vk3 + vk4)/6.0;
}

output_data::output_data(){
  init();
}

void output_data::init(){
  FILE *fp = fopen(FILE_PATH,"w");
    fprintf(fp, "time,");
  for(int i = 0; i < N_particle; i++) {
    fprintf(fp, "x%d,y%d,z%d,", i, i, i);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

void output_data::write(Particle_list const &P_list, double t){
  FILE *fp = fopen(FILE_PATH,"a");
  fprintf(fp, "%e,", t);
  for(int i = 0; i < N_particle; i++) {
    fprintf(fp, "%e,%e,%e,", P_list[i].r.x, P_list[i].r.y, P_list[i].r.z);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

void output_data::calc_time(int t){
  int day,hour,minute,second;
  day = t/(60*60*24);
  t = t%(60*60*24);
  hour = t/(60*60);
  t = t%(60*60);
  minute = t/60;
  t = t%60;
  second = t;
  printf("This calculatioin took %d day %d hour %d minute %d second\n\n", day, hour, minute, second);
}

