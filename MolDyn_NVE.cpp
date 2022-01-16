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
#include <ostream>
#include <fstream>      // Stream class to both read and write from/to files.
#include <iomanip>
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;
int main()
{
    Input(); //Inizialization
    int nconf = 1;
    ofstream Epot, Ekin, Temp, Etot, Pres;
    Epot.open("output_epot.dat");
    Ekin.open("output_ekin.dat");
    Temp.open("output_temp.dat");
    Etot.open("output_etot.dat");
    Pres.open("output_pres.dat");
    
    for (int iblk = 1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for (int istep = 1; istep <= (nstep / nblk); ++istep)
        {
            Move();
            Measure();
                Epot << stima_pot  << endl;
                Ekin << stima_kin  << endl;
                Temp << stima_temp << endl;
                Etot << stima_etot << endl;
                Pres << stima_pres << endl;
            Accumulate(); //Update block averages
            if (istep % 10 == 0) {
                //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
                nconf += 1;
            }
        }
        Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();
    return 0;
}

void Input(void) { //Prepare all stuff for the simulation
    ifstream ReadInput, ReadConf, ReadOldConf;
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
    vol = (double)npart / rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol, 1.0 / 3.0);
    cout << "Edge of the simulation box = " << box << endl;
  
    ReadInput >> rcut;

    //Correzioni di coda per energia potenziale e viriale
    vtail = (8.0 * pi * rho) / (9.0 * pow(rcut, 9)) - (8.0 * pi * rho) / (3.0 * pow(rcut, 3));
    ptail = (32.0 * pi * rho) / (9.0 * pow(rcut, 9)) - (16.0 * pi * rho) / (3.0 * pow(rcut, 3));
    cout << "Tail correction for the potential energy = " << vtail << endl;
    cout << "Tail correction for the virial           = " << ptail << endl;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> nblk;
    ReadInput >> iprint;


    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    ReadInput.close();

    //Prepare array for measurements
     iv = 0; //Potential energy
     ik = 1; //Kinetic energy
     ie = 2; //Total energy
     it = 3; //Temperature
     iw = 4; //Virial
     igofr = 5;
     nbins = 100;
     bin_size = (box / 2.0) / (double)nbins;
     n_props = 5; //Number of observables
     n_props = n_props + nbins;

    //Read initial configuration
    if (restart == 1) {
        cout << "Restart enabled" << endl;
        cout << "Read initial configuration from file old.0 " << endl;
        cout << "Read configuration at time t-dt from file old.final " << endl << endl;
        ReadConf.open("old.0");  //Reading configuration at actual time t
        ReadOldConf.open("old.final"); //Reading configuration at time t-dt
        for (int i = 0; i < npart; ++i) {
            ReadConf >> x[i] >> y[i] >> z[i];
            ReadOldConf >> xold[i] >> yold[i] >> zold[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;

            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        ReadConf.close();
        ReadOldConf.close();
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
    //Prepare random velocities in equilibration phase
    if (restart == 0) {
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = { 0.0, 0.0, 0.0 };
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

    if (restart == 1) {
        cout << "Prepare velocities with center of mass velocity equal to zero using Verlet's algorithm " << endl << endl;
        Move();
        double sumv2 = 0.0, fs;
        for (int i = 0; i < npart; ++i) {
            vx[i] = (x[i] - xold[i]) / 2 * delta;        //v(t+dt/2) = (x(t+dt) - x(t)) / 2dt
            vy[i] = (x[i] - xold[i]) / 2 * delta;
            vz[i] = (x[i] - xold[i]) / 2 * delta;

            sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        }
        sumv2 /= (double)npart;  
        fs = sqrt(3 * temp / sumv2);                              // fs = velocity scale factor 
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

void Move(void) { //Move particles with Verlet algorithm
    double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

    for (int i = 0; i < npart; ++i) { //Force acting on particle i
        fx[i] = Force(i, 0);
        fy[i] = Force(i, 1);
        fz[i] = Force(i, 2);
    }

    for (int i = 0; i < npart; ++i) { //Verlet integration scheme

        xnew = Pbc(2.0 * x[i] - xold[i] + fx[i] * pow(delta, 2));
        ynew = Pbc(2.0 * y[i] - yold[i] + fy[i] * pow(delta, 2));
        znew = Pbc(2.0 * z[i] - zold[i] + fz[i] * pow(delta, 2));

        vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
        vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
        vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);

        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];

        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
    }
    return;
}

double Force(int ip, int idir) { //Compute forces as -Grad_ip V(r)
    double f = 0.0;
    double dvec[3], dr;

    for (int i = 0; i < npart; ++i) {
        if (i != ip) {
            dvec[0] = Pbc(x[ip] - x[i]);  // distance ip-i in pbc
            dvec[1] = Pbc(y[ip] - y[i]);
            dvec[2] = Pbc(z[ip] - z[i]);

            dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
            dr = sqrt(dr);

            if (dr < rcut) {
                f += dvec[idir] * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8)); // -Grad_ip V(r)
            }
        }
    }

    return f;
}

void Measure() { //Properties measurement
    double v, t, w, vij, wij;
    double dx, dy, dz, dr;
    v = 0.0; //reset observables
    t = 0.0;
    w = 0.0;
    //reset the hystogram of g(r)
    for (int k = igofr; k < igofr + nbins; ++k) walker[k] = 0.0;
    //cycle over pairs of particles
    for (int i = 0; i < npart - 1; ++i) {
        for (int j = i + 1; j < npart; ++j) {

            if (restart == 0) {
                dx = Pbc(xold[i] - xold[j]); // here I use old configurations [old = r(t)]
                dy = Pbc(yold[i] - yold[j]); // to be compatible with EKin which uses v(t)
                dz = Pbc(zold[i] - zold[j]); // => EPot should be computed with r(t)
                dr = dx * dx + dy * dy + dz * dz;
                dr = sqrt(dr);
            }

            else if (restart == 1) {
                dx = Pbc(x[i] - x[j]); // here I use actual configurations
                dy = Pbc(y[i] - y[j]); // to be compatible with EKin which uses v(t)
                dz = Pbc(z[i] - z[j]); 
                dr = dx * dx + dy * dy + dz * dz;
                dr = sqrt(dr);
                for (int k = igofr; k < igofr + nbins; ++k) {
                    if (dr >= (k - 5) * bin_size && dr < (k - 4) * bin_size) {        //Parto da k-5 perché gofr parte da 5 e voglio anche valutare i primi bin
                        walker[k] += 2;
                    }
                }
            }

            if (dr < rcut) {
                vij = 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
                wij = 1.0 / pow(dr, 12) - 0.5 / pow(dr, 6);

                //Potential energy
                v += vij;
                w += wij;
            }
        }
    }
    walker[iv] = 4.0 * v;
    walker[iw] = 48.0 * w / 3.0;
    //Kinetic energy
    for (int i = 0; i < npart; ++i) t += (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

    walker[ik] = 0.5 * t;
    walker[it] = (2.0 / 3.0) * walker[ik];
    walker[ie] = walker[ik] + walker[iv];

    //stima_pot = walker[iv] / (double)npart;
    stima_pot = walker[iv] / (double)npart + vtail; //Potential energy per particle with tail correction
    stima_kin = walker[ik] / (double)npart; //Kinetic energy per particle
    stima_temp = walker[it] / (double)npart; //Temperature
    stima_etot = walker[ie] / (double)npart; //Total energy per particle
    //stima_pres = rho * temp + walker[iw] / vol;
    stima_pres = rho * temp + (walker[iw] + npart * ptail) / vol; //Total pressure with tail correction

    return;
}

void Reset(int iblk) //Reset block averages
{

    if (iblk == 1)
    {
        for (int i = 0; i < n_props; ++i)
        {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for (int i = 0; i < n_props; ++i)
    {
        blk_av[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}


void Accumulate(void) //Update block averages
{

    for (int i = 0; i < n_props; ++i)
    {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

    double r;
    ofstream Gofr, Gave, Ave_Epot, Ave_Ekin, Ave_Etot, Ave_Temp, Ave_Pres;
    const int wd = 12;

    cout << "Block number " << iblk << endl;
    
    Ave_Epot.open("ave_epot.dat", ios::app);
    Ave_Ekin.open("ave_ekin.dat", ios::app);
    Ave_Etot.open("ave_etot.dat", ios::app);
    Ave_Temp.open("ave_temp.dat", ios::app);
    Ave_Pres.open("ave_pres.dat", ios::app);
    Gofr.open("output.gofr.0", ios::app);
    Gave.open("output.gave.0", ios::app);
 

    stima_pot = blk_av[iv] / blk_norm / (double)npart + vtail; //Potential energy with tail correction
    //stima_pot = blk_av[iv] / blk_norm / (double)npart;    //Without tail correction
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot * stima_pot;
    err_pot = Error(glob_av[iv], glob_av2[iv], iblk);

    stima_kin = blk_av[ik] / blk_norm / (double)npart; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin * stima_kin;
    err_kin = Error(glob_av[ik], glob_av2[ik], iblk);

    stima_temp = blk_av[it] / blk_norm / (double)npart; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp * stima_temp;
    err_temp = Error(glob_av[it], glob_av2[it], iblk);

    stima_etot = blk_av[ie] / blk_norm / (double)npart; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot * stima_etot;
    err_etot = Error(glob_av[ie], glob_av2[ie], iblk);

    stima_pres = rho * temp + (blk_av[iw] / blk_norm + ptail * (double)npart) / vol; //Pressure with tail correction
    //stima_pres = rho * temp + (blk_av[iw] / blk_norm) / vol;  //Without tail correction
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres * stima_pres;
    err_press = Error(glob_av[iw], glob_av2[iw], iblk);

    for (int i = 5; i < n_props; i++) {
        delta_V = pow((i - 4) * bin_size, 3) - pow((i - 5) * bin_size, 3);
        stima_g = blk_av[i] / blk_norm / (rho * npart * (4 * M_PI / 3) * delta_V);
        glob_av[i] += stima_g;
        glob_av2[i] += stima_g * stima_g;
        err_gdir = Error(glob_av[i], glob_av2[i], iblk);
        r = (i - 5) * bin_size + bin_size / 2;
        Gofr << setw(wd) << iblk << setw(wd) << r << "\t " << glob_av[i] / (double)iblk << endl;
        if(iblk == nblk) Gave << setw(wd) << r << "\t " << glob_av[i] / (double)iblk << "\t " << err_gdir << endl;
    }

    Ave_Epot << setw(wd) << iblk << setw(wd) << stima_pot << setw(wd) << glob_av[iv] / (double)iblk << "\t" << err_pot << endl;
    Ave_Pres << setw(wd) << iblk << setw(wd) << stima_pres << setw(wd) << glob_av[iw] / (double)iblk << "\t" << err_press << endl;
    Ave_Ekin << setw(wd) << iblk << setw(wd) << stima_kin << setw(wd) << glob_av[ik] / (double)iblk << "\t" << err_kin << endl;
    Ave_Temp << setw(wd) << iblk << setw(wd) << stima_temp << setw(wd) << glob_av[it] / (double)iblk << "\t" << err_temp << endl;
    Ave_Etot << setw(wd) << iblk << setw(wd) << stima_etot << setw(wd) << glob_av[ie] / (double)iblk << "\t" << err_etot << endl;
    cout << "----------------------------" << endl << endl;

    Ave_Epot.close();
    Ave_Ekin.close();
    Ave_Temp.close();
    Ave_Etot.close();
    Ave_Pres.close();
    Gofr.close();
    Gave.close();
}

void ConfFinal(void) { //Write final configuration
    ofstream WriteConf;
    ofstream WriteOut;
    WriteOut.open("old.0"); //Writing config at time t
    WriteConf.open("old.final"); //Writing configuration at time t-dt

    cout << "Print previous configuration to old.final " << endl;
    cout << "Print final configuration to old.0" << endl << endl;

    for (int i = 0; i < npart; ++i) {
        WriteOut << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
        WriteConf << xold[i] / box << "   " << yold[i] / box << "   " << zold[i] / box << endl;
    }
    WriteConf.close();
    WriteOut.close();
    return;
}


void ConfXYZ(int nconf) { //Write configuration in .xyz format
    ofstream WriteXYZ;

    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for (int i = 0; i < npart; ++i) {
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " << Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
    }
    WriteXYZ.close();
}

double Pbc(double r) {  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r / box);
}

double Error(double sum, double sum2, int iblk)
{
    if (iblk == 1) return 0.0;
    else return sqrt((sum2 / (double)iblk - pow(sum / (double)iblk, 2)) / (double)(iblk - 1));
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

