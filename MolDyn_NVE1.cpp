/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#define _USE_MATH_DEFINES 
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  int N = 100;
  int L = nstep / N;
  int n = 0;
  double r = 0;
  double stima_g = 0;
  double delta_V = 0;
  sumprog_ekin = 0;
  sumprog_epot = 0;
  sumprog_etot = 0;
  sumprog_T = 0;
  mean_ekin = 0;
  mean_epot = 0;
  mean_etot = 0;
  mean_T = 0;
  cum_mean_ekin = 0;
  cum_mean_epot = 0;
  cum_mean_etot = 0;
  mean2_ekin = 0;
  mean2_epot = 0;
  mean2_etot = 0;
  mean2_T = 0;
  cum_mean2_ekin = 0;
  cum_mean2_epot = 0;
  cum_mean2_etot = 0;
  std_ekin = 0;
  std_epot = 0;
  std_etot = 0;
  std_T = 0;
  igofr = 2;
  nbins = 100;
  bin_size = (box / 2.0) / (double)nbins;

  ofstream Epot, Ekin, Etot, Temp, ave_ekin, ave_epot, ave_etot, ave_temp, Gave, Gofr;
  Epot.open("output_epot.dat");
  Ekin.open("output_ekin.dat");
  Temp.open("output_temp.dat");
  Etot.open("output_etot.dat");
  ave_ekin.open("ave_ekin.dat");
  ave_epot.open("ave_epot.dat");
  ave_etot.open("ave_etot.dat");
  ave_temp.open("ave_temp.dat");
  Gofr.open("output.gofr.0");
  Gave.open("output.gave.0");

  for(int j = 1; j <= N; ++j){
      sumprog_ekin = 0;
      sumprog_epot = 0;
      sumprog_etot = 0;
      sumprog_T = 0;
      n = 0;      
      for (int istep = 1; istep <= L; ++istep) {
         // if (istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
          Move();           //Move particles with Verlet algorithm
          if (istep%10 == 0) {
              Measure();     //Properties measurement
              Epot << stima_pot << endl;
              Ekin << stima_kin << endl;
              Temp << stima_temp << endl;
              Etot << stima_etot << endl;
              sumprog_ekin += stima_kin;            //Calcolo delle somme progressive dei valori di energia e T istantanei per particella
              sumprog_epot += stima_pot;
              sumprog_etot += stima_etot;
              sumprog_T += stima_temp;
              for (int i = 0; i < n_props; ++i) {
                  sumprog_g[i] += walker[i];
              }
              //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
              nconf += 1;
              n++;
          }
      }
      mean_ekin += sumprog_ekin / n;        //Calcolo della somma delle medie e delle medie quadratiche
      mean_epot += sumprog_epot / n;
      mean_etot += sumprog_etot / n;
      mean_T += sumprog_T / n;
      mean2_ekin += pow(sumprog_ekin / n, 2);
      mean2_epot += pow(sumprog_epot / n, 2);
      mean2_etot += pow(sumprog_etot / n, 2);
      mean2_T += pow(sumprog_T / n, 2);

      cum_mean_ekin = mean_ekin / j;        //Calcolo delle medie cumulative e delle medie quadratiche cumulative
      cum_mean_epot = mean_epot / j;
      cum_mean_etot = mean_etot / j;
      cum_mean_T = mean_T / j;
      cum_mean2_ekin = mean2_ekin / j;
      cum_mean2_epot = mean2_epot / j;
      cum_mean2_etot = mean2_etot / j;
      cum_mean2_T = mean2_T / j;

      for (int i = 0; i < n_props; i++) {
          delta_V = pow((i + 1) * bin_size, 3) - pow(i * bin_size, 3);
          stima_g = sumprog_g[i] / (n * rho * npart * (4 * M_PI / 3) * delta_V);
          cum_mean_g[i] += stima_g;
          cum_mean2_g[i] += (stima_g * stima_g);
          cum_mean_g[i] /= (double)j;
          cum_mean2_g[i] /= (double)j;
         
          if (j - 1 == 0) {
              std_g[i] = 0;
          }
          else {
              std_g[i] = sqrt(((cum_mean2_g[i]) - pow(cum_mean_g[i], 2)) / (j - 1));
          }         
          r = i * bin_size + bin_size / 2;
          Gofr << j << " " << r << "\t " << stima_g << endl;
          Gave << j << " " << r << "\t " << cum_mean_g[i]  << "\t " << std_g[i] << endl;

      }
      if (j == 1)
      {
          for (int i = 0; i < n_props; ++i)
          {
              cum_mean_g[i] = 0;
              cum_mean2_g[i] = 0;
          }
      }

      for (int i = 0; i < n_props; ++i)
      {
          sumprog_g[i] = 0;
          std_g[i] = 0;
      }

      if (j - 1 == 0) {                     //Calcolo della deviazione standard per la valutazione dell'errore sulle misure
          std_ekin = 0;
          std_epot = 0;
          std_etot = 0;
          std_T = 0;
      }
      else {
          std_ekin = sqrt((cum_mean2_ekin - pow(cum_mean_ekin, 2)) / (j - 1));
          std_epot = sqrt((cum_mean2_epot - pow(cum_mean_epot, 2)) / (j - 1));
          std_etot = sqrt((cum_mean2_etot - pow(cum_mean_etot, 2)) / (j - 1));
          std_T = sqrt((cum_mean2_T - pow(cum_mean_T, 2)) / (j - 1));
      }
       ave_ekin << j << " " << cum_mean_ekin << " " << std_ekin << endl;
       ave_epot << j << " " << cum_mean_epot << " " << std_epot << endl;
       ave_etot << j << " " << cum_mean_etot << " " << std_etot << endl;
       ave_temp << j << " " << cum_mean_T << " " << std_T << endl;     
  }

  ConfFinal();   //Write final configuration to restart Epot.close();
  Ekin.close();
  Epot.close();
  Temp.close();
  Etot.close();
  ave_ekin.close();
  ave_epot.close();
  ave_etot.close();
  ave_temp.close();
  Gave.close();
  Gofr.close();

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf, ReadOldConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

 ReadInput >> restart;
  if (restart == 1) {
      cout << "Restart enabled" << endl;
  }

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
 

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
 // iv = 0; //Potential energy
 // ik = 1; //Kinetic energy
 // ie = 2; //Total energy
 // it = 3; //Temperature
 // n_props = 4; //Number of observables

//Read initial configuration
  if (restart == 1) {
      cout << "Read previous time's configuration from file config.final " << endl << endl;
      ReadOldConf.open("config.final");
      for (int i = 0; i < npart; ++i) {
          ReadOldConf >> xold[i] >> yold[i] >> zold[i];
          xold[i] = xold[i] * box;
          yold[i] = yold[i] * box;
          zold[i] = zold[i] * box;
      }
      ReadOldConf.close();
      cout << "Read initial configuration from file config.final " << endl << endl;
      ReadConf.open("config.final");
      for (int i = 0; i < npart; ++i) {
          ReadConf >> x[i] >> y[i] >> z[i];
          x[i] = x[i] * box;
          y[i] = y[i] * box;
          z[i] = z[i] * box;
      }
      ReadConf.close();
  }
  else {
      cout << "Restart disabled: No old configuration available" << endl << endl;
      cout << "Read initial configuration from file config.0 " << endl << endl;
      ReadConf.open("config.0");
      for (int i = 0; i < npart; ++i) {
          ReadConf >> x[i] >> y[i] >> z[i];
          x[i] = x[i] * box;
          y[i] = y[i] * box;
          z[i] = z[i] * box;
      }
      ReadConf.close();
  }
  

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   if (restart == 0) {
       for (int i = 0; i < npart; ++i) {
           vx[i] = rand() / double(RAND_MAX) - 0.5;
           vy[i] = rand() / double(RAND_MAX) - 0.5;
           vz[i] = rand() / double(RAND_MAX) - 0.5;

           sumv[0] += vx[i];
           sumv[1] += vy[i];
           sumv[2] += vz[i];
       }
       for (int idim = 0; idim < 3; ++idim) sumv[idim] /= (double)npart;
       double sumv2 = 0.0, fs;
       for (int i = 0; i < npart; ++i) {
           vx[i] = vx[i] - sumv[0];
           vy[i] = vy[i] - sumv[1];
           vz[i] = vz[i] - sumv[2];

           sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
       }
       sumv2 /= (double)npart;

       fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 

       for (int i = 0; i < npart; ++i) {
           vx[i] *= fs;
           vy[i] *= fs;
           vz[i] *= fs;

           xold[i] = Pbc(x[i] - vx[i] * delta);
           yold[i] = Pbc(y[i] - vy[i] * delta);
           zold[i] = Pbc(z[i] - vz[i] * delta);
       }

   }
     
   if(restart == 1) {
       Move();
       double sumv2 = 0.0, fs;
       for (int i = 0; i < npart; ++i) {
           vx[i] = (x[i] - xold[i]) / 2 * delta;        //v(t+dt/2) = (x(t+dt) - x(t)) / 2dt
           vy[i] = (x[i] - xold[i]) / 2 * delta;
           vz[i] = (x[i] - xold[i]) / 2 * delta;

           sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
       }
       sumv2 /= (double)npart;
       //double T_new = sumv2 * 1. / 3;  
       fs = sqrt(3 * temp/ sumv2);                              // fs = velocity scale factor 
       for (int i = 0; i < npart; ++i) {
           vx[i] *= fs;
           vy[i] *= fs;
           vz[i] *= fs;          
           xold[i] = Pbc(x[i] - vx[i] * delta);
           yold[i] = Pbc(y[i] - vy[i] * delta);
           zold[i] = Pbc(z[i] - vz[i] * delta);
       }
   } 
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
    
    //if (restart == 0) {
        vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
        vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
        vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);
   // }
   
    /*else if(restart == 1) {
       
        vx_old[i] = Pbc(x[i] - xold[i]) /  delta;
        vy_old[i] = Pbc(y[i] - yold[i]) / delta;
        vz_old[i] = Pbc(z[i] - zold[i]) / delta;

        vx[i] = Pbc(vx_old[i] + delta * fx[i]);
        vy[i] = Pbc(vy_old[i] + delta * fy[i]);
        vz[i] = Pbc(vz_old[i] + delta * fz[i]);
       vx[i] = Pbc((xnew - x[i]) / delta);
        vy[i] = Pbc((ynew - y[i]) / delta);
        vz[i] = Pbc((znew - z[i]) / delta);

    }*/

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  //int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
 /* ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);*/

  v = 0.0; //reset observables
  t = 0.0;
  //reset the hystogram of g(r)
  for (int k = 0; k < nbins; ++k) walker[k] = 0.0;
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

        if (restart == 0) {
            dx = Pbc(xold[i] - xold[j]); // here I use old configurations [old = r(t)]
            dy = Pbc(yold[i] - yold[j]); // to be compatible with EKin which uses v(t)
            dz = Pbc(zold[i] - zold[j]); // => EPot should be computed with r(t)
            dr = dx * dx + dy * dy + dz * dz;
            dr = sqrt(dr);            
        }
        
        else if (restart == 1) {
            dx = Pbc(x[i] - x[j]); // here I use old configurations [old = r(t)]
            dy = Pbc(y[i] - y[j]); // to be compatible with EKin which uses v(t)
            dz = Pbc(z[i] - z[j]); // => EPot should be computed with r(t)
            dr = dx * dx + dy * dy + dz * dz;
            dr = sqrt(dr); 
            for (int k = 0; k < nbins; ++k) {
                if (dr >= k * bin_size && dr < (k + 1) * bin_size) {        //Parto da k-2 perché gofr parte da 2 e voglio anche valutare i primi bin
                    walker[k] += 2;
                }
            }
        }

        if (dr < rcut) {
            vij = 4.0 / pow(dr, 12) - 4.0 / pow(dr, 6);

            //Potential energy
            v += vij;
        }
       
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    /*Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();*/

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
  ofstream WriteOldConf;
  cout << "Print final configuration to file config.final " << endl << endl;
  cout << "Print final old configuration to file old.final " << endl << endl;
  WriteConf.open("config.final");
  WriteOldConf.open("old.final");
      for (int i = 0; i < npart; ++i) {
          WriteConf << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
          WriteOldConf << xold[i] / box << "   " << yold[i] / box << "   " << zold[i] / box << endl;
      }
      WriteConf.close();
      WriteOldConf.close();  
      return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
