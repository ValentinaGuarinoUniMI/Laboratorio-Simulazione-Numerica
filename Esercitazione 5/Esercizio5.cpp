#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm> 
#include "random.h"

using namespace std;

//Esercizio 5: Campionamento delle densità di probabilità delle funzioni d'onda per l'atomo di idrogeno relative allo stato 1s e allo 
// stato eccitato 2p. Per campionare le posizioni dalle suddette densità viene utilizzato l'algoritmo di Metropolis, utilizzando come
// trial transition probability T(x|y) sia una probabilità uniforme sia una probabilità gaussiana. Le posizioni campionate vengono usate per stimare
// il valore d'aspettazione del raggio nei due diversi stati, il quale verrà poi confrontato con il valore d'aspettazione trovato analiticamente.

int main() {
	//Inizializzazione delle variabili
	Random randgen;
	int seed[4];
	int p1, p2;
	int M = 1000000;
	int N = 200; 
	int L = M / N;
	int dim = 2;
	double delta = 1.2;
	double a0 = 1;		//Il raggio di Bohr viene impostato ad 1
	double A = 0;
	double A_precision = 0;
	double r = 0;
	double sumprog_r = 0;
	double mean_r = 0;
	double mean2_r = 0;
	double cumavr_r = 0;
	double cumavr2_r = 0;
	double std = 0;
	double x[dim] = { 0 };			//Nuove posizioni
	double y[dim] = { 1 };		//Vecchie posizioni: y = 1		//Partenza lontana dall'origine: y=500
	double p_old1s = 0;
	double p_new1s = 0;
	double p_old2p = 0;
	double p_new2p = 0;
	double sigma = 0.75;
	ofstream meanr1s, meanr1s_gauss, meanr2p, meanr2p_gauss;
	ofstream config1s_unif, config1s_gauss, config2p_unif, config2p_gauss;

	ifstream Primes("Primes");
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

	meanr1s.open("meanr1s_unif.dat");
	meanr2p.open("meanr2p_unif.dat");
	meanr1s_gauss.open("meanr1s_gauss.dat");
	meanr2p_gauss.open("meanr2p_gauss.dat");
	config1s_unif.open("config1s_unif.dat");
	config1s_gauss.open("config1s_gauss.dat");
	config2p_unif.open("config2p_unif.dat");
	config2p_gauss.open("config2p_gauss.dat");

//Stato 1s: campionamento con T(x|y) uniforme
	for (int i = 1; i <= N; i++) {				//Blocking method per eliminare le correlazioni date dall'algoritmo di Metropolis
		sumprog_r = 0;
		for (int j = 0; j < L; j++) {
			x[0] = randgen.Rannyu(y[0] - delta, y[0] + delta);
			x[1] = randgen.Rannyu(y[1] - delta, y[1] + delta);
			x[2] = randgen.Rannyu(y[2] - delta, y[2] + delta);
//Densità di probabilità per lo stato 1s			
			p_old1s = (pow(a0, -3) / M_PI) * exp((-2 * sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2])) / a0);
			p_new1s = (pow(a0, -3) / M_PI) * exp((-2 * sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])) / a0);
			
			A = min(1., p_new1s / p_old1s);			//Applicazione del criterio di accettazione
			if (A == 1) {
				y[0] = x[0];
				y[1] = x[1];
				y[2] = x[2];
				A_precision++;
			}
			else {
				r = randgen.Rannyu();
				if (A >= r) {
					y[0] = x[0];
					y[1] = x[1];
					y[2] = x[2];
					A_precision++;
				}
			}
			sumprog_r += sqrt(pow(y[0],2) + pow(y[1],2) + pow(y[2],2));
			//config1s_unif << x[0] << " " << x[1] << " " << x[2] << endl;  //Togliere il commento per avere la stampa delle posizioni campionate
		}
		mean_r += sumprog_r / L;
		mean2_r += pow(sumprog_r / L, 2);
		cumavr_r = mean_r / i;
		cumavr2_r = mean2_r / i;
		
		if ((i - 1) == 0) {
			std = 0;
		}
		else {
			std = sqrt((cumavr2_r - pow(cumavr_r, 2)) / (i - 1));
		}

		meanr1s << i << " " << cumavr_r << " " << std << endl;
		
	}
	A_precision /= M;
	cout << "Precision for state 1s:" << " " << A_precision << "\n";	//Calcolo della precisione dell'algoritmo: ci si aspetta un'accettazione di circa 0.5
	meanr1s.close();
	config1s_unif.close();

	//Azzeramento delle variabili 
	x[dim] = { 0 };
	y[dim] = { 1 };
	sumprog_r = 0;
	mean_r = 0;
	mean2_r = 0;
	cumavr_r = 0;
	cumavr2_r = 0;
	std = 0;
	//Stato 1s: campionamento con T(x|y) gaussiana multivariata
	for (int i = 1; i <= N; i++) {
		sumprog_r = 0;
		for (int j = 0; j < L; j++) {
			x[0] = randgen.Gauss(y[0], sigma);
			x[1] = randgen.Gauss(y[1], sigma);
			x[2] = randgen.Gauss(y[2], sigma);

			p_old1s = (pow(a0, -3) / M_PI) * exp((-2 * sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2])) / a0);
			p_new1s = (pow(a0, -3) / M_PI) * exp((-2 * sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])) / a0);

			A = min(1., p_new1s / p_old1s);
			if (A == 1) {
				y[0] = x[0];
				y[1] = x[1];
				y[2] = x[2];
				A_precision++;
			}
			else {
				r = randgen.Rannyu();
				if (A >= r) {
					y[0] = x[0];
					y[1] = x[1];
					y[2] = x[2];
					A_precision++;
				}
			}
			sumprog_r += sqrt(pow(y[0], 2) + pow(y[1], 2) + pow(y[2], 2));
			//config1s_gauss << x[0] << " " << x[1] << " " << x[2] << endl;    //Togliere il commento per avere la stampa delle posizioni campionate
		}
		mean_r += sumprog_r / L;
		mean2_r += pow(sumprog_r / L, 2);
		cumavr_r = mean_r / i;
		cumavr2_r = mean2_r / i;

		if ((i - 1) == 0) {
			std = 0;
		}
		else {
			std = sqrt((cumavr2_r - pow(cumavr_r, 2)) / (i - 1));
		}

		meanr1s_gauss << i << " " << cumavr_r << " " << std << endl;
	}
	A_precision /= M;
	cout << "Precision for state 1s with gaussian T:" << " " << A_precision << "\n";
	meanr1s_gauss.close();
	config1s_gauss.close();

//Azzeramento variabili
	x[dim] = { 0 };
	y[dim] = { 1 };
	sumprog_r = 0;
	mean_r = 0;
	mean2_r = 0;
	cumavr_r = 0;
	cumavr2_r = 0;
	std = 0;
	delta = 3;
	a0 = 1;
//Stato 2p: campionamento con T(x|y) uniforme
	for (int i = 1; i <= N; i++) {
		sumprog_r = 0;
		for (int j = 0; j < L; j++) {
			x[0] = randgen.Rannyu(y[0] - delta, y[0] + delta);
			x[1] = randgen.Rannyu(y[1] - delta, y[1] + delta);
			x[2] = randgen.Rannyu(y[2] - delta, y[2] + delta);

//Densità di probabilità per lo stato 2p
			p_old2p = (pow(a0, -5) / (32 * M_PI)) * pow(y[2], 2) * exp((-sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2])) / a0);
			p_new2p = (pow(a0, -5) / (32 * M_PI)) * pow(x[2], 2) * exp((-sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])) / a0);

			A = min(1., p_new2p / p_old2p);

			if (A == 1) {
				y[0] = x[0];
				y[1] = x[1];
				y[2] = x[2];
				A_precision++;
			}
			else {
				r = randgen.Rannyu();
				if (A >= r) {
					y[0] = x[0];
					y[1] = x[1];
					y[2] = x[2];
					A_precision++;
				}
			}
			sumprog_r += sqrt(pow(y[0], 2) + pow(y[1], 2) + pow(y[2], 2));
			//config2p_unif << x[0] << " " << x[1] << " " << x[2] << endl;	//Togliere il commento per avere la stampa delle posizioni campionate
		}
		mean_r += sumprog_r / L;
		mean2_r += pow(sumprog_r / L, 2);
		cumavr_r = mean_r / i;
		cumavr2_r = mean2_r / i;

		if ((i - 1) == 0) {
			std = 0;
		}
		else {
			std = sqrt((cumavr2_r - pow(cumavr_r, 2)) / (i - 1));
		}
		meanr2p << i << " " << cumavr_r << " " << std << endl;
		
	}
	A_precision /= M;
	cout << "Precision for state 2p:" << " " << A_precision << "\n";
	meanr2p.close();
	config2p_unif.close();
//Azzeramento variabili
	x[dim] = { 0 };
	y[dim] = { 1};
	sumprog_r = 0;
	mean_r = 0;
	mean2_r = 0;
	cumavr_r = 0;
	cumavr2_r = 0;
	std = 0;
	sigma = 1.8;
//Stato 2p: campionamento con T(x|y) gaussiana multivariata
	for (int i = 1; i <= N; i++) {
		sumprog_r = 0;
		for (int j = 0; j < L; j++) {
			x[0] = randgen.Gauss(y[0], sigma);
			x[1] = randgen.Gauss(y[1], sigma);
			x[2] = randgen.Gauss(y[2], sigma);

			p_old2p = (pow(a0, -5) / (32 * M_PI)) * pow(y[2], 2) * exp((-sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2])) / a0);
			p_new2p = (pow(a0, -5) / (32 * M_PI)) * pow(x[2], 2) * exp((-sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])) / a0);

			A = min(1., p_new2p / p_old2p);

			if (A == 1) {
				y[0] = x[0];
				y[1] = x[1];
				y[2] = x[2];
				A_precision++;
			}
			else {
				r = randgen.Rannyu();
				if (A >= r) {
					y[0] = x[0];
					y[1] = x[1];
					y[2] = x[2];
					A_precision++;
				}
			}
			sumprog_r += sqrt(pow(y[0], 2) + pow(y[1], 2) + pow(y[2], 2));
			//config2p_gauss << x[0] << " " << x[1] << " " << x[2] << endl;		//Togliere il commento per avere la stampa delle posizioni campionate
		}
		mean_r += sumprog_r / L;
		mean2_r += pow(sumprog_r / L, 2);
		cumavr_r = mean_r / i;
		cumavr2_r = mean2_r / i;

		if ((i - 1) == 0) {
			std = 0;
		}
		else {
			std = sqrt((cumavr2_r - pow(cumavr_r, 2)) / (i - 1));
		}
		meanr2p_gauss << i << " " << cumavr_r << " " << std << endl;
		
	}
	A_precision /= M;
	cout << "Precision for state 2p with gaussian T:" << " " << A_precision;
	meanr2p_gauss.close();
	config2p_gauss.close();
	return 0;
}