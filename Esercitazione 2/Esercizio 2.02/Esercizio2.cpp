#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "random.h"

using namespace std;
/*
* Esercizio 2.02: implementazione di un random walk 3D in un reticolo cubico nel caso discreto e nel caso continuo. Simulo un cammino con 100 passi, 
* in ognuno dei quali il walker può andare avanti o indietro, e faccio 10^4 ripetizioni del processo. Analizzo l'andamento della radice della 
* media dello spostamento quadratico rispetto al numero di passi; l'errore sulla quantità viene calcolato tramite propagazione degli errori.
*/

int main() {                            //Dichiarazione e inizializzazione delle variabili
    Random randgen;
    int seed[4];
    int p1, p2;
    double r0, w;
    double theta, phi;
    double x1 = 0;
    double y1 = 0;
    double z1 = 0;
    double x = 0;
    double y = 0;
    double z = 0;
    int throws = 10000;
    int steps = 100;
    int a = 1;
    double rho=0;
    double r =  0;
    double meanr[steps] = { 0 };
    double meanr2[steps] = { 0 };
    double walk[steps] = { 0 };
    double std[steps] = { 0 };
    double meanrho[steps] = { 0 };
    double meanrho2[steps] = { 0 };
    double walk_cont[steps] = { 0 };
    double std_cont[steps] = { 0 };
    double mean_meanr[steps] = { 0 };
    double mean_meanrho[steps] = { 0 };
    double mean_meanr2[steps] = { 0 };
    double mean_meanrho2[steps] = { 0 };
    ofstream rwdiscr;
    ofstream rwcont;

    ifstream Primes("Primes");                              //Comandi per la scelta del seme utilizzato nel generatore dei numeri pseudorandom
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    }
    else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                randgen.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }
    else cerr << "PROBLEM: Unable to open seed.in" << endl;

    rwdiscr.open("rwdiscr.dat");
    rwcont.open("rwcont.dat");

    for (int i = 0; i < throws; i++) {               //Ciclo sulle ripetizioni
        x = 0;                                       //Inizializzo a zero le coordinate
        y = 0;
        z = 0;
        x1 = 0;
        y1 = 0;
        z1 = 0;
        for (int j = 0; j < steps-1; j++) {           //Ciclo sul numero dei passi          
                r0 = randgen.Rannyu();                //Scelta della direzione: il RW può avvenire in avanti o indietro
                if (r0 < 0.5) {
                    r0 = -1;
                }
                else {
                    r0 = 1;
                }
                w = randgen.Rannyu();                 //Caso discreto: la scelta della direzione del passo avviene tramite campionamento casuale
                theta = randgen.UnifAngle();          //Caso continuo: il valore delle coordinate si trova campionando gli angoli theta e phi
                phi = randgen.Rannyu(0, 2*M_PI);
                if (w < 0.33) {                      //Passi in una direzione nel caso discreto
                    x +=  a * r0;
                }
                else if (w > 0.33 && w < 0.66) {
                    y +=  a * r0 ;                   
                }
                else {
                    z +=  a * r0 ;                   
                }   
                x1 += a * (sin(theta) * cos(phi));   //Passi nel caso continuo
                y1 += a * (sin(theta) * sin(phi));
                z1 += a * cos(theta);

                rho = x1 * x1 + y1 * y1 + z1 * z1;  //Calcolo delle distanze al quadrato: r=r^2, rho=rho^2
                r = x * x + y * y + z * z;
                meanr[j + 1] += r;                  //Poiché all'istante iniziale ci si trova nell'origine, i vettori delle medie vengono riempiti
                meanrho[j + 1] += rho;              //dal secondo passo
                meanr2[j + 1] += r * r;
                meanrho2[j + 1] += rho * rho;                
        }     
    }    
    for(int i = 0; i < steps; i++) {                //Ciclo sul numero di passi: per avere le medie divido per il numero totale di ripetizioni
        meanr[i] = meanr[i] / throws;
        meanrho[i] = meanrho[i] / throws;
        meanr2[i] = meanr2[i] / throws;
        meanrho2[i] = meanrho2[i] / throws;
        mean_meanr[i] = meanr[i] / throws;
        mean_meanrho[i] = meanrho[i] / throws;
        mean_meanr2[i] = meanr2[i] / throws;
        mean_meanrho2[i] = meanrho2[i] / throws;
        walk[i] = sqrt(meanr[i]);
        walk_cont[i] = sqrt(meanrho[i]);
        if (i == 0) {                               //Calcolo dell'errore utilizzando la propagazione degli errori
            std[i] = 0;
            std_cont[i] = 0;
        }
        else {
            std[i] = 1./2 * pow(meanr[i], -1./2) * sqrt((mean_meanr2[i] - mean_meanr[i] * mean_meanr[i]) / i);
            std_cont[i] = 1./2 * pow(meanrho[i], -1./2) * sqrt((mean_meanrho2[i] - mean_meanrho[i] * mean_meanrho[i]) / i);
        }
        rwdiscr << i << " " << walk[i] << " " << std[i] << endl;
        rwcont << i << " " << walk_cont[i] << " " << std_cont[i] << endl;
    }        
    rwcont.close();
    rwdiscr.close();
    randgen.SaveSeed();
    return 0;
}

