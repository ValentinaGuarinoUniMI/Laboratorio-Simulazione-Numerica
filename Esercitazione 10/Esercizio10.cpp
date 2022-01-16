#define _USE_MATH_DEFINES
#include<iostream>
#include<stdio.h>
#include<cmath>
#include<random>
#include<vector>
#include<algorithm>
#include<ctime>
#include<cstdlib>
#include<fstream>
#include "SimAnn.h"

using namespace std;

int main() {

	int length_path = 32;	//Numero di città
	double scaling_factor_high = 0.5;			//T > 2
	double scaling_factor_low = 0.9999;			//T < 2
	double T = 30;
	double A = 0;
	double old_dist, old, new_dist, old_dist_square, new_dist_square = 0; 
	double beta = 1. / T;
	double r;
	int n;
	int n_iter_high = 10;
	int n_iter_low = 50;
	vector<vector<double>> position;
	vector<vector<double>> position_square;
	random_device rand;
	default_random_engine gen(rand());
	uniform_real_distribution<> unif_distr(0.0, 1.0);
	
	
	ofstream output_old;
	ofstream output_best;
	ofstream output_L;
	ofstream output_meanL;
	ofstream output_old_square;
	ofstream output_best_square;
	ofstream output_L_square;
	ofstream output_meanL_square;
	
	output_old.open("oldpath.dat");
	output_best.open("bestpath.dat");
	output_L.open("L1.dat");
	output_old_square.open("oldpath_square.dat");
	output_best_square.open("bestpath_square.dat");
	output_L_square.open("L1_square.dat");

//Ricerca percorso migliore per città poste su una circonferenza	
	SetPosition_Circumference(position, length_path);

	Path my_path(position, length_path);
	Path new_path(position, length_path);
	
	cout << "Initial path and distance on circumference " << endl;
	my_path.PrintPath();
	my_path.PrintDistance();
	old_dist = my_path.GetDistance();
	old = old_dist;

	for(int i = 0; i < length_path; i++) {
		output_old << my_path.GetCity()[i].GetX() << " " << my_path.GetCity()[i].GetY() << " ";
		output_old << my_path.GetCity()[i].GetAllele() << " " << my_path.GetDistance() << endl;
		
	}
	output_old << my_path.GetCity()[0].GetX() << " " << my_path.GetCity()[0].GetY() << " ";
	output_old << my_path.GetCity()[0].GetAllele() << " " << my_path.GetDistance() << endl;
	
	n = 0;
	cout << "Initial T = " << T << endl;
	output_L << n << " " << old_dist << endl;
	
	do {
		beta = 1. / T;
		
		for(int i = 0; i < n_iter_high; i++) {
			
			new_path.Pair_Permutation();	//Applico le varie mosse al percorso e valuto di volta in volta l'accettazione
			new_dist = new_path.GetDistance();
			A = min(1., exp(-beta*(new_dist - old_dist)));
			if(A == 1) {
				old_dist = new_dist;
				my_path = new_path;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist = new_dist;
					my_path = new_path;
				}
				
			}
			
			
			new_path = my_path;
			
			new_path.Inversion_Permutation();	
			new_dist = new_path.GetDistance();
			A = min(1., exp(-beta*(new_dist - old_dist)));
			if(A == 1) {
				old_dist = new_dist;
				my_path = new_path;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist = new_dist;
					my_path = new_path;
				}
				
			}
			
			new_path = my_path;
			
			new_path.Contiguous_Permutation();	
			new_dist = new_path.GetDistance();
			A = min(1., exp(-beta*(new_dist - old_dist)));
			if(A == 1) {
				old_dist = new_dist;
				my_path = new_path;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist = new_dist;
					my_path = new_path;
				}
				
			}
			
			new_path = my_path;
			
			new_path.Shift_Mutation();	
			new_dist = new_path.GetDistance();
			A = min(1., exp(-beta*(new_dist - old_dist)));
			if(A == 1) {
				old_dist = new_dist;
				my_path = new_path;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist = new_dist;
					my_path = new_path;
				}
				
			}
			
			new_path = my_path;	
				
		}
	
		T = T * scaling_factor_high;	
		n++;	
		if(n % 1000 == 0) {
			output_L << n << " " << old_dist << endl;
			cout << "Iteration " << n << endl;
		}
		
	} while(T > 2);
	
	do {
		beta = 1. / T;
		
		for(int i = 0; i < n_iter_low; i++) {
			
			new_path.Pair_Permutation();	//Applico le varie mosse al percorso e valuto di volta in volta l'accettazione
			new_dist = new_path.GetDistance();
			A = min(1., exp(-beta*(new_dist - old_dist)));
			if(A == 1) {
				old_dist = new_dist;
				my_path = new_path;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist = new_dist;
					my_path = new_path;
				}
				
			}
			
			new_path = my_path;
			
			new_path.Inversion_Permutation();	
			new_dist = new_path.GetDistance();
			A = min(1., exp(-beta*(new_dist - old_dist)));
			if(A == 1) {
				old_dist = new_dist;
				my_path = new_path;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist = new_dist;
					my_path = new_path;
				}
				
			}
			
			new_path = my_path;
			
			new_path.Contiguous_Permutation();	
			new_dist = new_path.GetDistance();
			A = min(1., exp(-beta*(new_dist - old_dist)));
			if(A == 1) {
				old_dist = new_dist;
				my_path = new_path;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist = new_dist;
					my_path = new_path;
				}
				
			}
			
			new_path = my_path;
			
			new_path.Shift_Mutation();	
			new_dist = new_path.GetDistance();
			A = min(1., exp(-beta*(new_dist - old_dist)));
			if(A == 1) {
				old_dist = new_dist;
				my_path = new_path;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist = new_dist;
					my_path = new_path;
				}
				
			}
			
			new_path = my_path;	
					
				
		}
	
		T = T * scaling_factor_low;	
		n++;	
		if(n % 1000 == 0) {
			output_L << n << " " << old_dist << endl;
			cout << "Iteration " << n << endl;
		}
		
	} while(T > 0.005 && T < 2);
	
	
	cout << "Final best path and distance on circumference " << endl;
	
	my_path.PrintPath();
	my_path.PrintDistance();
	
	cout << "Final T = " << T << endl;
	
	for(int i = 0; i < length_path; i++) {
		output_best << my_path.GetCity()[i].GetX() << " " << my_path.GetCity()[i].GetY() << " ";
		output_best << my_path.GetCity()[i].GetAllele() << " " << my_path.GetDistance() << " ";
		output_best << old << endl;
		
	}
	output_best << my_path.GetCity()[0].GetX() << " " << my_path.GetCity()[0].GetY() << " ";	//Condizione di ritorno alla città iniziale
	output_best << my_path.GetCity()[0].GetAllele() << " " << my_path.GetDistance() << " ";
	output_best << old << endl;
	
	output_L.close();
	output_best.close();
	output_old.close();
	
//Ricerca percorso migliore per città poste all'interno di un quadrato
	SetPosition_Square(position_square, length_path);
	
	Path my_path_square(position_square, length_path);
	Path new_path_square(position, length_path);	
		
	cout << "Initial best path and distance on square " << endl;
	my_path_square.PrintPath();
	my_path_square.PrintDistance();
	
	
	old_dist_square = my_path_square.GetDistance();	
	old = old_dist_square;
	
	for(int i = 0; i < length_path; i++) {
		output_old_square << my_path_square.GetCity()[i].GetX() << " " << my_path_square.GetCity()[i].GetY() << " ";
		output_old_square << my_path_square.GetCity()[i].GetAllele() << " " << my_path_square.GetDistance() << endl;
		
	}
	output_old_square << my_path_square.GetCity()[0].GetX() << " " << my_path_square.GetCity()[0].GetY() << " ";
	output_old_square << my_path_square.GetCity()[0].GetAllele() << " " << my_path_square.GetDistance() << endl;
	
	
	T = 30;
	n = 0;
	cout << "Initial T = " << T << endl;
	output_L_square << n << " " << old_dist_square << endl;
	
	do {
	
		beta = 1. / T;
		for(int i = 0; i < n_iter_high; i++) {
	
			new_path_square.Pair_Permutation();	//Applico le varie mosse al percorso e valuto di volta in volta l'accettazione
			new_dist_square = new_path_square.GetDistance();
			
			A = min(1., exp(-beta*(new_dist_square - old_dist_square)));
			if(A == 1) {
				old_dist_square = new_dist_square;
				my_path_square = new_path_square;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist_square = new_dist_square;
					my_path_square = new_path_square;
				}
				
			}
			
			new_path_square = my_path_square;
			
			new_path_square.Inversion_Permutation();	
			new_dist_square = new_path_square.GetDistance();
			
			A = min(1., exp(-beta*(new_dist_square - old_dist_square)));
			if(A == 1) {
				old_dist_square = new_dist_square;
				my_path_square = new_path_square;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist_square = new_dist_square;
					my_path_square = new_path_square;
				}
				
			}
			
			new_path_square = my_path_square;
			
			new_path_square.Contiguous_Permutation();	
			new_dist_square = new_path_square.GetDistance();
			
			A = min(1., exp(-beta*(new_dist_square - old_dist_square)));
			if(A == 1) {
				old_dist_square = new_dist_square;
				my_path_square = new_path_square;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist_square = new_dist_square;
					my_path_square = new_path_square;
				}
				
			}
			
			new_path_square = my_path_square;
			
			new_path_square.Shift_Mutation();	
			new_dist_square = new_path_square.GetDistance();
			
			A = min(1., exp(-beta*(new_dist_square - old_dist_square)));
			if(A == 1) {
				old_dist_square = new_dist_square;
				my_path_square = new_path_square;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist_square = new_dist_square;
					my_path_square = new_path_square;
				}
				
			}
			
			new_path_square = my_path_square;
			
		}
		
		T = scaling_factor_high * T;
		n++;
		if(n % 1000 == 0) {
			output_L_square << n << " " << old_dist_square << endl;
			cout << "Iteration " << n << endl;
		}
		
	} while(T > 2);
	
	do {
	
		beta = 1. / T;
		for(int i = 0; i < n_iter_low; i++) {
	
			new_path_square.Pair_Permutation();	//Applico le varie mosse al percorso e valuto di volta in volta l'accettazione
			new_dist_square = new_path_square.GetDistance();
			
			A = min(1., exp(-beta*(new_dist_square - old_dist_square)));
			if(A == 1) {
				old_dist_square = new_dist_square;
				my_path_square = new_path_square;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist_square = new_dist_square;
					my_path_square = new_path_square;
				}
				
			}
			
			new_path_square = my_path_square;
			
			new_path_square.Inversion_Permutation();	
			new_dist_square = new_path_square.GetDistance();
			
			A = min(1., exp(-beta*(new_dist_square - old_dist_square)));
			if(A == 1) {
				old_dist_square = new_dist_square;
				my_path_square = new_path_square;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist_square = new_dist_square;
					my_path_square = new_path_square;
				}
				
			}
			
			new_path_square = my_path_square;
			
			new_path_square.Contiguous_Permutation();	
			new_dist_square = new_path_square.GetDistance();
			
			A = min(1., exp(-beta*(new_dist_square - old_dist_square)));
			if(A == 1) {
				old_dist_square = new_dist_square;
				my_path_square = new_path_square;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist_square = new_dist_square;
					my_path_square = new_path_square;
				}
				
			}
			
			new_path_square = my_path_square;
			
			new_path_square.Shift_Mutation();	
			new_dist_square = new_path_square.GetDistance();
			
			A = min(1., exp(-beta*(new_dist_square - old_dist_square)));
			if(A == 1) {
				old_dist_square = new_dist_square;
				my_path_square = new_path_square;
			}
			else {
				r = unif_distr(gen);
				if(A > r) {
					old_dist_square = new_dist_square;
					my_path_square = new_path_square;
				}
				
			}
			
			new_path_square = my_path_square;
			
		}
		
		T = scaling_factor_low * T;
		n++;
		if(n % 1000 == 0) {
			output_L_square << n << " " << old_dist_square << endl;
			cout << "Iteration " << n << endl;
		}
		
	} while(T > 0.005 && T < 2);
	
		
	
	cout << "Final best path and distance on square " << endl;
	
	my_path_square.PrintPath();
	my_path_square.PrintDistance();
	
	cout << "Final T = " << T << endl;
	
	for(int i = 0; i < length_path; i++) {
		output_best_square << my_path_square.GetCity()[i].GetX() << " " << my_path_square.GetCity()[i].GetY() << " ";
		output_best_square << my_path_square.GetCity()[i].GetAllele() << " " << my_path_square.GetDistance() << " ";
		output_best_square << old<< endl;
		
	}
	output_best_square << my_path_square.GetCity()[0].GetX() << " " << my_path_square.GetCity()[0].GetY() << " ";	//Condizione di ritorno alla città iniziale
	output_best_square << my_path_square.GetCity()[0].GetAllele() << " " << my_path_square.GetDistance() << " ";
	output_best_square << old << endl;
	
	output_L_square.close();
	output_best_square.close();
	output_old_square.close();
	
	return 0;
}
