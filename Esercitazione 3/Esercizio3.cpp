#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "random.h"

using namespace std;
/*Esercizio 3: calcolo di un European call-option price, C[S(0),0], e put-option price, P[S(0),0], tramite sampling Monte-Carlo
di un Moto Browniano Geometrico. La call e la put vengono valutate sia campionando in modo diretto l'asset finale del prezzo
S(T) = S(0) * exp((r - 0.5 * vol^2) * T + vol * W(T), dove W(T)~N(0,T) è un processo di Wiener, sia campionando il GBM(r,vol^2) tramite discretizzazione
dell'asset dei prezzi, S(t_i+1) = S(t_i) * exp((r - 0.5 * vol^2) * delta_t + vol * Z_i * sqrt(delta_t)) dove delta_t = t_i+1 - t_i e Z_i~N(0,1).
Utilizzando il metodo del data blocking per il calcolo delle incertezze statistiche, si confrontano i risultati ottenuti al variare del numero di blocchi
con i risultati analitici di Black-Scholes per un Plain Vanilla Option Pricing. */

int main(){

	Random randgen;				//Inizializzazione delle variabili
	int M = 10000;
	int N = 100;
	int L = M / N;
	int seed[4];
	int p1, p2;
	
	int S_0 = 100;
	double K = 100;
	int T = 1;
	double r = 0.1;
	double vol = 0.25;
	int n_steps = 100;
	double S_T;
	double S_t = S_0;
	double C;
	double P;
	double z;
	double w;
	double sum_prog_call = 0;
	double sum_prog_put = 0;
	double delta_t = T * 1. / n_steps;			//Per il sampling dell'asset del prezzo discretizzato divido l'intervallo [0,T] in 100 sottointervalli
	double mean_call = 0;
	double mean_put = 0;
	double mean2_put = 0;
	double mean2_call =  0;
	double cumavr_call = 0;
	double cumavr2_call = 0;
	double cumavr_put = 0;
	double cumavr2_put = 0;
	double std_call = 0;
	double std_put = 0;

	ofstream call_dir;
	ofstream put_dir;
	ofstream call_discr;
	ofstream put_discr;

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

	call_dir.open("call_dir.dat");
	put_dir.open("put_dir.dat");
	call_discr.open("call_discr.dat");
	put_discr.open("put_discr.dat");

	//Sampling nel caso di campionamento diretto dell'asset price finale al tempo T utilizzando data blocking
	for (int i = 1; i <= N; i++) {
		sum_prog_call = 0;
		sum_prog_put = 0;
		for (int j = 0; j < L; j++) {
			z = randgen.Gauss(0, 1);
			S_T = S_0 * exp((r - 0.5 * vol * vol) * T + vol * z * sqrt(T));
			C = exp(-r * T) * max(S_T - K, 0.0);
			P = exp(-r * T) * max(K - S_T, 0.0);
			sum_prog_call += C;
			sum_prog_put += P;
		}
		mean_call += sum_prog_call / L;
		mean2_call += pow(sum_prog_call / L, 2);
		mean_put += sum_prog_put / L;
		mean2_put += pow(sum_prog_put / L, 2);

		cumavr_call = mean_call / i;
		cumavr2_call = mean2_call / i;
		cumavr_put = mean_put / i;
		cumavr2_put = mean2_put / i;
		if (i - 1 == 0) {
			std_call = 0;
			std_put = 0;
		}
		else {
			std_call = sqrt((cumavr2_call - pow(cumavr_call, 2)) / (i - 1));
			std_put = sqrt((cumavr2_put - pow(cumavr_put, 2)) / (i - 1));
		}
		call_dir << i << " " << cumavr_call << " " << std_call << endl;
		put_dir << i << " " << cumavr_put << " " << std_put << endl;
	}
	C = 0;						//Rinizializzazione delle variabili
	P = 0;
	mean_call = 0;
    mean_put = 0;
	mean2_put = 0;
	mean2_call = 0;
	cumavr_call = 0;
	cumavr2_call = 0;
	cumavr_put = 0;
	cumavr2_put = 0;
	std_call = 0;
	std_put = 0;

	//Sampling per un asset price S(t) discretizzato utilizzando il metodo di data blocking
	for (int i = 1; i <= N; i++) {
		sum_prog_call = 0;
		sum_prog_put = 0;
		for (int j = 0; j < L; j++) {
			S_t = S_0;
			for (int t = 1; t < N; t++) {
				w = randgen.Gauss(0, 1);
				S_t = S_t * exp((r - 0.5 * vol * vol) * delta_t + vol * w * sqrt(delta_t));		
			}			
			C = exp(-r * T) * max(0.0, S_t -K);
			P = exp(-r * T) * max(K - S_t, 0.0);
			sum_prog_call += C;
			sum_prog_put += P;
		}
		mean_call += sum_prog_call / L;
		mean2_call += pow(sum_prog_call / L, 2);
		mean_put += sum_prog_put / L;
		mean2_put += pow(sum_prog_put / L, 2);

		cumavr_call = mean_call / i;
		cumavr2_call = mean2_call / i;
		cumavr_put = mean_put / i;
		cumavr2_put = mean2_put / i;
		if (i - 1 == 0) {
			std_call = 0;
			std_put = 0;
		}
		else {
			std_call = sqrt((cumavr2_call - pow(cumavr_call, 2)) / (i-1));
			std_put = sqrt((cumavr2_put - pow(cumavr_put, 2)) / (i - 1));
		}
		call_discr << i << " " << cumavr_call << " " << std_call << endl;
		put_discr << i << " " << cumavr_put << " " << std_put << endl;
	}
	call_dir.close();
	put_dir.close();
	call_discr.close();
	put_discr.close();
	return 0;
}
