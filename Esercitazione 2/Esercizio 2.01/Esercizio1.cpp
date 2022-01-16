#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "random.h"
/*Esercizio 2.01: calcolo di un integrale utilizzando metodi Monte Carlo. Viene valutata la differenza tra il valore vero dell'integrale e quello ottenuto
* con la tecnica Monte Carlo (come media di una funzione G opportunamente valutata) in due casi diversi: campionamento uniforme tra [0,1] e importance sampling. 
* Come funzione nell'importance sampling viene testata sia f(x) = 2(1 - x) sia f(x) = 3/2*sqrt(1 - x). Si utilizza il metodo a blocchi.
*/
using namespace std;

int main() {				//Dichiarazione e inizializzazoine variabili utilizzate

	Random randgen;
	int seed[4];
	int p1, p2;
	int M = 100000;
	int N = 100;
	int L = M / N;
	double r;
	double x;
	//double z;
	double sum_G = 0;
	double sum_G_is = 0;
	double mean_G = 0;
	double mean2_G = 0;
	double mean_G_is = 0;
	double mean2_G_is = 0;
	double cumavr_G = 0;
	double cumavr2_G = 0;
	double cumavr_G_is = 0;
	double cumavr2_G_is = 0;
	double std = 0;
	double std_is = 0;

	ofstream unif;
	ofstream impsampl;

	ifstream Primes("Primes");		//Comandi per la scelta del seme utilizzato nel generatore dei numeri pseudorandom
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

	unif.open("unif.dat");
	impsampl.open("impsampl.dat");

	for (int i = 1; i <= N; i++) {					//Ciclo sui blocchi
		sum_G = 0;									//Inizializzazione a zero delle variabili somma per il ciclo successivo
		sum_G_is = 0;
		for (int j = 0; j < L; j++) {				//Ciclo per ogni esperimento nel blocco
			r = randgen.Rannyu();
			x = randgen.ImpSampl();					//Numeri random distribuiti secondo le f(x) di importance sampling scelte. Ottenuti tramite 
			//z = randgen.ImpSampl2();				//inversione della funzione cumulativa nel random.cpp
			sum_G += (M_PI/2)*cos((M_PI * r) / 2);	//Caso uniforme
			sum_G_is += float((M_PI/3) * cos((M_PI * x) / 2) / sqrt(1 - x));	//Caso importance sampling 1:  f(x) = 3/2*sqrt(1 - x)
			//sum_G_is += float((M_PI / 4) * cos((M_PI * z) / 2) / (1 - z));    //Caso importance sampling 2:  f(x) = 2(1 - x)
		}
		mean_G += sum_G / L;					//Calcolo delle medie
		mean_G_is += sum_G_is / L;
		mean2_G += pow(sum_G / L, 2);
		mean2_G_is += pow(sum_G_is / L, 2);
		cumavr_G = mean_G / i;
		cumavr2_G = mean2_G / i;
		cumavr_G_is = mean_G_is / i;
		cumavr2_G_is = mean2_G_is / i;
		if (i - 1 == 0) {								//Calcolo la deviazione standard per le medie e le varianze cumulate di ogni blocco
			std = 0;
			std_is = 0;
		}
		else {
			std = sqrt((cumavr2_G - pow(cumavr_G, 2)) / (i - 1));
			std_is = sqrt((cumavr2_G_is - pow(cumavr_G_is, 2)) / (i - 1));
		}
		unif << i << " " <<cumavr_G << " " << std << endl;					//STampa su file dei valori ottenuti
		impsampl << i << " " << cumavr_G_is << " " << std_is << endl;
	}
	impsampl.close();
	unif.close();
	randgen.SaveSeed();
	return 0;
}