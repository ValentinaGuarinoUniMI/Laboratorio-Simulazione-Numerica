#include <iostream>
#include <cmath>
#include <fstream>
#include "random.h"

using namespace std;
/*
Esercizio 1.02: verifica del CLT attraverso la generazione di RV uniformi,esponenziali e lorenziane. Nel generatore di numeri casuali si sono implementate
le funzioni per il campionamento di numeri pseudorandom che seguono distribuzioni esponenziali e lorenziane grazie al metodo dell'inversione della
funzione cumulativa. Vengono generati 10^4 numeri casuali e viene calcolata la loro somma S_N = (1/N) * sum(x_i) per diversi valori di N={1,2,10,100}.
Si valuta come al crescere di N le distribuzioni uniforme ed esponenziale convergano ad una gaussiana, mentre la distribuzione lorenziana converga ad una
lorenziana.
*/
int main() {

	int p1, p2;					//Dichiarazione e inizializzazione variabili
	int lambda = 1;
	int gamma = 1;
	int mean = 0;
	int N;
	int throws = 10000;
	double unif;
	double lor;
	double exp;
	double sum_unif = 0;
	double sum_exp = 0;
	double sum_lor = 0;

	Random randgen;

	int seed[4];

	ofstream unif_out;
	ofstream lor_out;
	ofstream exp_out;

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

	unif_out.open("unif_out.dat");
	exp_out.open("exp_out.dat");
	lor_out.open("lor_out.dat");

	for (N = 0; N < 4; N++) {							//Ciclo sui diversi possibili valori della variabile N={1,2,10,100}
		if (N == 0) {
			for (int i = 0; i < throws; i++) {			//N=1: genero 10^4 numeri casuali dalle diverse distribuzioni. In questo caso S_N = x_i.
				unif = randgen.Rannyu();
				exp = randgen.Exp(lambda);
				lor = randgen.Lorentz(gamma, mean);

				unif_out << unif << " ";
				exp_out << exp << " ";
				lor_out << lor << " ";
			}
			unif_out << endl;
			exp_out << endl;
			lor_out << endl;
		}		
		else if (N == 1) {							//N=2: genero 2 RV e ne calcolo la somma S_2 = 1/2(x_1+x_2); lo faccio per 10^4 volte
			for (int i = 0; i < throws; i++) {
				sum_unif = 0;
				sum_exp = 0;
				sum_lor = 0;
				for(int j=0; j<2; j++) {			//Ciclo per sommare le N=2 RV
					unif = randgen.Rannyu();
					exp = randgen.Exp(lambda);
					lor = randgen.Lorentz(gamma, mean);

					sum_unif += unif;
					sum_exp += exp;
					sum_lor += lor;
				}
				sum_unif /= 2;
				sum_exp /= 2;
				sum_lor /= 2;
				unif_out << sum_unif << " ";
				exp_out << sum_exp << " ";
				lor_out << sum_lor << " ";
			}
			unif_out << endl;
			exp_out << endl;
			lor_out << endl;
		}			
		else if (N == 2) {						//N=10: genero 10 RV e calcolo la variabile somma S_N; lo faccio per 10^4 volte
			for (int i = 0; i < throws; i++) {
				sum_unif = 0;
				sum_exp = 0;
				sum_lor = 0;
				for (int j = 0; j < 10; j++) {
					unif = randgen.Rannyu();
					exp = randgen.Exp(lambda);
					lor = randgen.Lorentz(gamma, mean);

					sum_unif += unif;
					sum_exp += exp;
					sum_lor += lor;
				}
				sum_unif /= 10;
				sum_exp /= 10;
				sum_lor /= 10;
				unif_out << sum_unif << " ";
				exp_out << sum_exp << " ";
				lor_out << sum_lor << " ";
			}
			unif_out << endl;
			exp_out << endl;
			lor_out << endl;
		}		
	    else if (N == 3) {								//N=100: genero 100 RV e calcolo la variabile somma S_N; lo faccio per 10^4 volte
			for (int i = 0; i < throws; i++) {
				sum_unif = 0;
				sum_exp = 0;
				sum_lor = 0;
				for (int j = 0; j < 100; j++) {
					unif = randgen.Rannyu();
					exp = randgen.Exp(lambda);
					lor = randgen.Lorentz(gamma, mean);

					sum_unif += unif;
					sum_exp += exp;
					sum_lor += lor;
				}
				sum_unif /= 100;
				sum_exp /= 100;
				sum_lor /= 100;
				unif_out << sum_unif << " ";
				exp_out << sum_exp << " ";
				lor_out << sum_lor << " ";
			}
			unif_out << endl;
			exp_out << endl;
			lor_out << endl;
		}
			
	}	
	unif_out.close();
	exp_out.close();
	lor_out.close();
	randgen.SaveSeed();
	return 0;
}
