#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm> 
#include "random.h"

using namespace std;
//Esercizio 8.01-8.02: si utilizza il principio variazionale per la ricerca del valor medio dell'energia di ground state del sistema dato. 
//L'energia minima a cui arriva il sistema, che equivale al valore più vicino dell'energia di ground state, si trova ottimizzando i parametri dai quali
//dipende la funzione trial utilizzata per lo studio del sistema. L'ottimizzazione viene svolta creando una griglia di parametri e calcolando il valor medio
//dell'energia relativo ad ogni coppia di parametri. Per la funzione trial scelta, pari alla somma di due gaussiane, i parametri da ottimizzare sono il 
//valor medio mu e la deviazione standard sigma. Viene inoltre generato un file contenente i punti dello spazio visitati dall'algoritmo di Metropolis,
//in modo da poter creare successivamente un istogramma relativo alla distribuzione di probabilità campionata.
int main() {	

	//Inizializzazione variabili
	Random randgen;
	int seed[4];
	int p1, p2;
	int M = pow(10, 6);
	int N = pow(10, 2);
	int Z = pow(10, 4);
	int L = M / N;
	double x = 1;
	double x_new = 0;
	double x_new_opt = 0;
	double delta = 2.7;
	double delta_opt = 2.5;
	double A_precision = 0;
	double A_precision_opt = 0;
	int n_rows = 0;
	double h_bar = 1;
	double m = 1;
	double Energy = 0;
	double E_min = 0;
	double sumprog_H = 0, mean_H = 0, mean2_H = 0, cumavr_H = 0, cumavr2_H = 0, std = 0;
	double mu_old = 1;
	double sigma_old = 0.7;
	double mu_min = 0;
	double sigma_min = 0;
	double mu, sigma;
	double V = pow(x, 4) - 2.5 * pow(x, 2);
	double A, r, E_old, E_new, mu_new, sigma_new, p_old, p_new;
	double psi_trial, kin_psi, pot_psi;
	ofstream mean_ene, varparam, out_ene, out_min, histo_psi;


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

	mean_ene.open("mean_ene.dat");
	varparam.open("varparam.dat");
	out_ene.open("out_ene.dat");
	out_min.open("out_min.dat");
	histo_psi.open("histo_psi.dat");

//Creazione griglia di parametri
	for (int i = 0; i < 100; i++) {
		sigma_old = 0.7;
		mu_new = mu_old - 0.005;
		for (int j = 0; j < 100; j++) {
			sigma_new = sigma_old - 0.001;
			varparam << mu_new << " " << sigma_new << endl;
			sigma_old = sigma_new;
		}
		mu_old = mu_new;
	}

	ifstream ReadInput;
	ReadInput.open("varparam.dat");

//Ottimizzazione: Ricerca dei parametri che minimizzino l'energia
	while (!ReadInput.eof()) {
		n_rows++;
		ReadInput >> mu >> sigma;		//Lettura da file delle 
		E_old = 0;
//Algoritmo di Metropolis per l'ottimizzazione dei parametri 
		for (int i = 0; i < Z; i++) {
			
			x_new_opt = randgen.Rannyu(x + delta_opt, x - delta_opt);

			p_old = pow(exp(-pow((x - mu), 2) / (2 * pow(sigma, 2))) + exp(-pow((x + mu), 2) / (2 * pow(sigma, 2))), 2);
			p_new = pow(exp(-pow((x_new_opt - mu), 2) / (2 * pow(sigma, 2))) + exp(-pow((x_new_opt + mu), 2) / (2 * pow(sigma, 2))), 2);

			A = min(1., p_new / p_old);		//Calcolo dell'accettazione
			if (A == 1) {
				x = x_new_opt;
				A_precision_opt++;
			}
			else {
				r = randgen.Rannyu();
				if (A >= r) {
					x = x_new_opt;
					A_precision_opt++;
				}
			}
			psi_trial = exp(-pow((x - mu), 2) / (2 * pow(sigma, 2))) + exp(-pow((x + mu), 2) / (2 * pow(sigma, 2)));
			V = pow(x, 4) - 2.5 * pow(x, 2);
			kin_psi = -(h_bar / (2 * m)) * (1. / pow(sigma, 2)) * (exp(-pow((x - mu), 2) / (2 * pow(sigma, 2))) * ((pow(x - mu, 2) / pow(sigma, 2)) - 1) + exp(-pow((x + mu), 2) / (2 * pow(sigma, 2))) * ((pow(x + mu, 2) / pow(sigma, 2)) - 1));
			pot_psi = V * psi_trial;
			E_old += (kin_psi + pot_psi) / psi_trial;
		}
		Energy = E_old / Z*1.0;
		out_ene << Energy << " " << mu << " " << sigma << endl;
		
		//Confronto le energie per ricercare i parametri con i quali ho l'energia minima
			if (Energy < E_min) {
				E_min = Energy;
				sigma_min = sigma;
				mu_min = mu;
			}
			out_min << E_min << " " << mu_min << " " << sigma_min << endl;
	}

	//Metodo a blocchi: ricerca dei valori medi cumulativi per l'energia e creazione dell'istogramma
	for (int i = 1; i <= N; i++) {
		sumprog_H = 0;
		for (int j = 0; j < L; j++) {

			x_new = randgen.Rannyu(x + delta, x - delta);

			p_old = pow(exp(-pow((x - mu_min), 2) / (2 * pow(sigma_min, 2))) + exp(-pow((x + mu_min), 2) / (2 * pow(sigma_min, 2))), 2);
			p_new = pow(exp(-pow((x_new - mu_min), 2) / (2 * pow(sigma_min, 2))) + exp(-pow((x_new + mu_min), 2) / (2 * pow(sigma_min, 2))), 2);

			A = min(1., p_new / p_old);
			if (A == 1) {
				x = x_new;
				A_precision++;
				histo_psi << x_new << endl;
			}
			else {
				r = randgen.Rannyu();
				if (A >= r) {
					x = x_new;
					A_precision++;
					histo_psi << x_new << endl;
				}
			}
			psi_trial = exp(-pow((x - mu_min), 2) / (2 * pow(sigma_min, 2))) + exp(-pow((x + mu_min), 2) / (2 * pow(sigma_min, 2)));
			V = pow(x, 4) - 2.5 * pow(x, 2);
			kin_psi = -(h_bar / (2 * m)) * (1. / pow(sigma_min, 2)) * (exp(-pow((x - mu_min), 2) / (2 * pow(sigma_min, 2))) * ((pow(x - mu_min, 2) / pow(sigma_min, 2)) - 1) + exp(-pow((x + mu_min), 2) / (2 * pow(sigma_min, 2))) * ((pow(x + mu_min, 2) / pow(sigma_min, 2)) - 1));
			pot_psi = V * psi_trial;
			E_new = (kin_psi + pot_psi) / psi_trial;

			sumprog_H += E_new;
		}

		mean_H += sumprog_H / L;
		mean2_H += pow(sumprog_H / L, 2);
		cumavr_H = mean_H / i;
		cumavr2_H = mean2_H / i;


		if ((i - 1) == 0) {
			std = 0;
		}
		else {
			std = sqrt((cumavr2_H - pow(cumavr_H, 2)) / (i - 1));
		}

		cout << "Blocco numero " << i << endl;
		mean_ene << i << " " << cumavr_H << " " << std << endl;
		//histo_psi << x[0] << " " << x[1] << " " << x[2] << endl;
	}
	A_precision /= M;
	A_precision_opt /= (Z * n_rows);
	cout << "Acceptance rate optimization's Metropolis " << A_precision_opt * 100 << "% with delta_opt = " << delta_opt << endl;
	cout << "Acceptance rate blocking's Metropolis " << A_precision * 100 << "% with delta = " << delta << endl;

	mean_ene.close();
	out_ene.close();
	out_min.close();
	histo_psi.close();

	return 0;
}