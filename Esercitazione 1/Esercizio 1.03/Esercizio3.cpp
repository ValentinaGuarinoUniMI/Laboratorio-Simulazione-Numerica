#include <iostream>
#include <cmath>
#include <fstream>
#include "random.h"

using namespace std;
/*
Esercizio 1.03: stima del valore di pi tramite l'esperimento di Buffon. Nel generatore di numeri casuali si è implementata una funzione per il campionamento
con distribuzione uniforme di un angolo tra [0,pi/2] utilizzando il metodo del rigetto. Calcolo la probabilità che l'ago intersechi una sbarra; ciò
avviene quando la distanza x tra il centro dell'ago e la sbarra è minore o uguale alla metà della lunghezza totale dell'ago per il seno 
dell'angolo relativo alla sua inclinazione. Campionando randomicamente l'angolo theta e l'inclinazione x, valutiamo quando vi è intersezione; la 
probabilità è valutata dunque attraverso la frequenza come N_lancitotali/N_successi.
Utilizzando il metodo a blocchi si calcola il valor medio di pi su ognuno degli N blocchi, per poi valutarne la media cumulata al variare dei blocchi.
*/

int main() {

	Random randgen;		//Dichiarazione e inizializzazione variabili
	ofstream outpi;
	int seed[4];
	int p1, p2;
	int M = 100000;
	int N = 100;
	int L = M / N;
	double x;
	double theta;
	int N_hit = 0;
	double d = 2.55;				//Scelgo d>l ma non d>>l
	int l = 1;
	double pi = 0;
	double avr_pi = 0;
	double avr2_pi = 0;
	double cumavr_pi = 0;
	double cumavr2_pi = 0;
	double std = 0;

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

	outpi.open("outpi.dat");						//Apertura del file di scrittura

	for (int i = 1; i <= N; i++) {					//Ciclo sui diversi N blocchi
		N_hit = 0;
		pi = 0;
		for (int j = 0; j < L; j++) {				//Ciclo su ognuno degli L tentativi all'interno di ogni blocco
			x = randgen.Rannyu();
			theta = randgen.UnifAngle();
			if (x <= (1. / 2) * sin(theta)) {		//Se la condizione è vera, l'ago interseca la sbarra; incremento il numero di successi
				N_hit += 1;
			}
			if (N_hit != 0) {
				pi += (2 * l * j) / (N_hit * d);		//Sommo i valori di pi ottenuti per ogni r,theta casuali
			}
		}
		avr_pi += pi / L;				//Calcolo media e media al quadrato del valore di pi
		avr2_pi += pow(pi / L, 2);
		cumavr_pi = avr_pi / i;			//Calcolo media cumulata
		cumavr2_pi = avr2_pi / i;
		if (i - 1 == 0) {								//Calcolo l'errore sulla misurazione
			std = 0;
		}
		else {
			std = sqrt((cumavr2_pi - pow(cumavr_pi, 2)) / (i-1));
		}
		outpi << i << " " << cumavr_pi << " " << std << endl;		//Stampa dei valori
	}
	outpi.close();
	randgen.SaveSeed();
	return 0;
}