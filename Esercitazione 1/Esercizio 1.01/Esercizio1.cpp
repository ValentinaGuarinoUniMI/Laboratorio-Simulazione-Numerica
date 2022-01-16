#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
/*Esercizio 1.01: valutiamo la bonta' del generatore di numeri pseudorandom vedendo lo scarto che si ha tra la media e la varianza calcolate per 
* M numeri casuali (divisi con raggruppamento a blocchi) rispetto ai valori noti di media e varianza per distribuzione uniforme. 
* Si valuta la frazione di numeri casuali che cade in un certo sottointervallo di [0,1], e si calcola il chi^2 rispetto al valore d'aspettazione che si 
* avrebbe  per una distribuzione binomiale (E(x) = np).
*/
int main() {

	Random randgen;		//Definisco un metodo appartenente alla classe Random, che genera numeri pseudo-random
	int seed[4];		//Seme iniziale utilizzato
	int M = 100000;		//Definizione e inizializzazione delle variabili utilizzate
	int N = 100;
	int L = M / N;
	int subint = 100;
	int throws = 10000;
	int exper = 100;
	double chi2;
	double r;
	int p1, p2;
	
	double avr = 0;
	double avr2 = 0;
	double var = 0;
	double var2 = 0;
	double cumavr = 0;
	double cumavr2 = 0;
	double cumvar = 0;
	double cumvar2 = 0;
	double sum_prog = 0;
	double sumvar_prog = 0;
	double std = 0;
	double std_var = 0;
	double n[subint] = { 0 };
	double inf = 0, sup = 0; 
	
	ofstream outavr;				//Comandi per abilitare scrittura su file
	ofstream outvar; 
	ofstream outchi;

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
	
//Punti 1 e 2: calcolo della media e della varianza di numeri pseudorandom su una distribuzione uniforme utilizzando metodo raggruppamento a blocchi

	outavr.open("outavr.dat");			
	outvar.open("outvar.dat");

	for (int i = 1; i <= N; i++) {		//Ciclo sui diversi blocchi N
		sum_prog = 0;					//Inizializzazione variabili somma
		sumvar_prog = 0;
		for (int j = 0; j < L; j++) {		//Ciclo su tutti gli L numeri casuali generati per ogni blocco N. Vengono incrementate le somme e vengono
			r = randgen.Rannyu();			//riempiti i vettori contenenti le medie trovate per ogni blocco. Si fa lo stesso per la varianza.
			sum_prog += r;
			sumvar_prog += pow(r - 0.5,2);
		}
		avr += sum_prog / L;
		avr2 += pow(sum_prog / L, 2);
		var += sumvar_prog / L;
		var2 += pow(sumvar_prog / L, 2);
		
		cumavr = avr / i;		//Al denominatore vi è i+1 in quanto i va da 0 a N-1
		cumavr2 = avr2 / i;
		cumvar = var / i;
		cumvar2 = var2 / i;
		if (i - 1 == 0) {								//Calcolo la deviazione standard per le medie e le varianze cumulate di ogni blocco
			std = 0;
			std_var = 0;
		}
		else {
			std = sqrt((cumavr2 - pow(cumavr, 2)) / (i-1));
			std_var = sqrt((cumvar2 - pow(cumvar, 2)) / (i-1));
		}
		outavr << i << " " << cumavr << " " << std << endl;		//Stampa dei valori per l'intervallo complessivo M
		outvar << i << " " << cumvar << " " << std_var << endl;
	}
	outavr.close();
	outvar.close();

//Punto 3: calcolo del chi^2. Divido [0,1] in 100 sottointervalli e faccio 100 esperimenti; per ogni esperimento genero 10^4 numeri casuali e valuto quanti
//ne cadono in ogni intervallo. Successivamente computo il chi^2.

	outchi.open("outchi.dat");

	for (int i = 0; i < exper; i++) {			//Ciclo sul numero di esperimenti
		double Xsqr;
		for (int j = 0; j < throws; j++) {		//Ciclo per la generazione di 10^4 numeri casuali
			r = randgen.Rannyu();
			for (int k = 0; k < subint; k++) {	//Ciclo sui diversi sotto-intervalli
				inf = k * 1. / subint;
				sup = (k + 1) * 1. / subint;
				if (r >= inf && r < sup) {		//Se il numero generato appartiene al sotto-intervallo considerato viene incrementato il contatore
					n[k] += 1;
				}
			}			
		}
		for (int l = 0; l < subint; l++) {		//Calcolo del chi^2 
			Xsqr = n[l] - (throws / subint);
			chi2 += (Xsqr * Xsqr) / (throws / subint);
		}
		outchi << i << " " << chi2 << endl;		//Stampa dei valori
		Xsqr = 0;
		chi2 = 0;
		for (int h = 0; h < subint; h++) {		//Svuotamento del contatore, lo reinizializzo per il nuovo esperimento
			n[h] = 0;
		}
	}
	randgen.SaveSeed();
	return 0;
}