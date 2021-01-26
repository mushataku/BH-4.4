#include "header.h"

/******************config********************/
// 0:重力のみ 1:Lorentzt力のみ(半径依存無し) 2:Lorentzt力のみ(半径依存有) 3:重力+Lorentz
const int FORCE = 3;
// 0:Euler 1:Runge-Kutta 2:Leap Flog
const int METHOD = 1;
/******************config********************/

/******************計算条件********************/
double r_ring = 6.0;
// const double beta = 0.0;
double beta = M_PI/10;
double beta_deg = beta/M_PI*180.0;
// スピンパラメータ
double a_spin = 1.0;
int N_particle = 50;
const int nt = 100001;
const int dn = 1000;
// 時間の基準とする公転速度の半径
const double r_Tg = 10.0;
const double tmax = 2.0*M_PI*std::pow(r_Tg, 1.5)*100.2;
// double tmax = M_PI/(2.0*a_spin)*10;
// double tmax = M_PI/(2.0*a_spin)*r_ring*r_ring*r_ring*3;
const double dt = tmax/double(nt-1);
/******************計算条件********************/

int main(int argc, char *argv[]){
  if(argc > 1){
    r_ring = atof(argv[1]);
    beta_deg = atof(argv[2]);
    a_spin = atof(argv[3]);
    beta = beta_deg*M_PI/180.0;
    // cast により自動的に切り捨て
    N_particle = 50*(r_ring/10);
    printf("Number of Particles = %d\n", N_particle);
  }
  // 時間計測
  clock_t start_t, end_t;
  start_t = time(NULL);

  double t = 0.0;
  Particle_list P_list(N_particle);

  // 出力
  P_list.output(t);

  printf("****************CALICULATION START****************\n");
  printf("r = %.1f, beta = %.1f°, a_spin = %.1f\n", r_ring, beta_deg, a_spin);
  printf("****************CALICULATION START****************\n");

  for(int i = 0; i < nt; i++) {
    update(P_list, t, dt);
    if(i%dn == 0) P_list.output(t);
  }

  printf("******************NORMAL END******************\n");
  end_t = time(NULL);
  // printf("This calculatioin took %ld second \n\n", end_t - start_t);
  calc_time((int)(end_t-start_t));
  return 0;
}

/******************関数********************/

void Particle_list::Particle::init(double r0, double phi){
  using namespace std;// ここだけ std:: を省略するために
  double v0, v_beta, v_beta_tmp, v_phi, tmp;
  v0 = (2.0*a_spin*cos(beta)+sqrt(4.0*a_spin*a_spin*cos(beta)*cos(beta)+r0*r0*r0))/(r0*r0);

  r = {r0*cos(phi), r0*cos(beta)*sin(phi), r0*sin(beta)*sin(phi)};
  v = {-v0*sin(phi), v0*cos(beta)*cos(phi), v0*sin(beta)*cos(phi)};
}

void Particle_list::init_particles(){
  double phi, r0, v0, x, y, z;
  for(int i = 0; i < N_particle; i++) {
    phi = 2.0*M_PI*i/N_particle;
    P_list[i].init(r_ring, phi);
  }
}

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
    if(METHOD == 0) euler(P_list.P_list[i].r, P_list.P_list[i].v, dt);
    if(METHOD == 1) runge_kutta(P_list.P_list[i].r, P_list.P_list[i].v, dt);
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

void Particle_list::init_file(){
  char DIR_PATH[50];
  sprintf(DIR_PATH, "./output/a%.1fbeta%.1f", a_spin, beta_deg);
  mkdir(DIR_PATH, 0777);
  sprintf(FILE_PATH, "%s/a%.1fbeta%.1fr%.1f.csv", DIR_PATH, a_spin, beta_deg, r_ring);
  FILE *fp = fopen(FILE_PATH,"w");
  fprintf(fp, "time,");
  for(int i = 0; i < N_particle; i++) {
    fprintf(fp, "x%d,y%d,z%d,r%d,E%d,vx%d,vy%d,vz%d,", i, i, i, i, i, i, i, i);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

void Particle_list::output(double t){
  FILE *fp = fopen(FILE_PATH, "a");
  fprintf(fp, "%e,", t);
  for(int i = 0; i < N_particle; i++) {
    fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,", P_list[i].r.x, P_list[i].r.y, 
      P_list[i].r.z, P_list[i].r.length(), P_list[i].Energy(),
      P_list[i].v.x, P_list[i].v.y, P_list[i].v.z);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

void calc_time(int t){
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

