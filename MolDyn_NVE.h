/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __NVE__
#define __NVE__

//parameters, observables
const int m_props=1000;
double vtail, ptail;
int n_props, iv, ik, it, ie, iw;
double sd;
double walker[m_props];
double delta_V;

// averages
double acc, att;
double blk_av[m_props], blk_norm;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres, stima_g;
double err_pot, err_kin, err_temp, err_etot, err_gdir, err_press;
//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part], vy[m_part], vz[m_part];
double vx_old[m_part], vy_old[m_part], vz_old[m_part];

// thermodynamical state
int npart,attempted,accepted;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, iprint, seed;
double delta;
double restart;

//pigreco
const double pi = 3.1415927;

//measurement of g(r)
double igofr;
double nbins;
double bin_size;
//functions
void Input(void);
void Move(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
#endif
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
