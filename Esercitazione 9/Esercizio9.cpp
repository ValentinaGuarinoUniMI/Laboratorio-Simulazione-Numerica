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
#include "GenAlg.h"

using namespace std;

int main() {

	int N_pop = 100;	//Grandezza della popolazione
	int length_path = 32;	//Numero di città
	int gen = 500;		//Numero di generazioni
	
	int best = 0;
	double dist, old_dist, old_dist_square = 0; 
	
	vector<vector<double>> position;
	vector<vector<double>> position_square;
	vector <Path> path;
	vector <Path> path_square;
	
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
	output_meanL.open("meanL1.dat");
	output_old_square.open("oldpath_square.dat");
	output_best_square.open("bestpath_square.dat");
	output_L_square.open("L1_square.dat");
	output_meanL_square.open("meanL1_square.dat");

//Ricerca percorso migliore per città poste su una circonferenza	

	//SetPosition_Circumference(position, length_path);
	SetPosition_FromFile(position, length_path);	//Valuto l'andamento al variare delle generazioni del percorso migliore iniziale di 								  una data simulazione
		
	for (int i = 0; i < N_pop; i++) {
		//RandomPath(position, length_path);	//Commentato per prendere come percorso iniziale uno già dato
		path.push_back(Path(position, length_path));
	}

	Fitness(path);
		
	/*for (auto i : path) {
		i.PrintPath();
		cout << "\n" << endl;
	      	i.PrintDistance();
	}*/
	cout << "Initial best path and distance on circumference " << endl;
	path[best].PrintPath();
	path[best].PrintDistance();
	old_dist = path[best].GetDistance();
	
	for(int i = 0; i < length_path; i++) {
		output_old << path[best].GetCity()[i].GetX() << " " << path[best].GetCity()[i].GetY() << " ";
		output_old << path[best].GetCity()[i].GetAllele() << " " << path[best].GetDistance() << endl;
	}
	output_old << path[best].GetCity()[0].GetX() << " " << path[best].GetCity()[0].GetY() << " ";
	output_old << path[best].GetCity()[0].GetAllele() << " " << path[best].GetDistance() << endl;
	
	NewGeneration newgen(path);
	
	for(int i = 0; i < gen; i++) {
		newgen.Evolution(N_pop);
		newgen.Mutation();
	
		output_L << i << " " << newgen.GetPath()[0].GetDistance() << endl;
		for(int j = 0; j < N_pop / 2; j++) {
			dist += newgen.GetPath()[j].GetDistance();
		}
		output_meanL << i << " " << dist * 1. / (1.*N_pop/2.) << endl;
		dist = 0;
	}
	
	cout << "Final best path and distance on circumference " << endl;
	newgen.Fitness();
	newgen.GetPath()[best].PrintPath();
	newgen.GetPath()[best].PrintDistance();
	for(int i = 0; i < length_path; i++) {
		output_best << newgen.GetPath()[best].GetCity()[i].GetX() << " " << newgen.GetPath()[best].GetCity()[i].GetY() << " ";
		output_best << newgen.GetPath()[best].GetCity()[i].GetAllele() << " " << newgen.GetPath()[best].GetDistance() << " ";
		output_best << old_dist << endl;
		
	}
	output_best << newgen.GetPath()[best].GetCity()[0].GetX() << " " << newgen.GetPath()[best].GetCity()[0].GetY() << " ";	//Condizione di ritorno alla città iniziale
	output_best << newgen.GetPath()[best].GetCity()[0].GetAllele() << " " << newgen.GetPath()[best].GetDistance() << " ";
	output_best << old_dist << endl;
	
	output_L.close();
	output_meanL.close();
	output_best.close();
	output_old.close();
	
//Ricerca percorso migliore per città poste all'interno di un quadrato

	//SetPosition_Square(position_square, length_path);
	SetPosition_FromFile_Square(position_square, length_path);	//Valuto l'andamento al variare delle generazioni del percorso 										  migliore iniziale di una data simulazione
	dist = 0;
	
	for (int i = 0; i < N_pop; i++) {
		//RandomPath(position_square, length_path);	//Commentato per prendere come percorso iniziale uno già dato
		path_square.push_back(Path(position_square, length_path));
	}

	Fitness(path_square);
		
	/*for (auto i : path) {
		i.PrintPath();
		cout << "\n" << endl;
	      	i.PrintDistance();
	}*/
	cout << "Initial best path and distance on square " << endl;
	path_square[best].PrintPath();
	path_square[best].PrintDistance();
	
	old_dist_square = path_square[best].GetDistance();
	
	for(int i = 0; i < length_path; i++) {
		output_old_square << path_square[best].GetCity()[i].GetX() << " " << path_square[best].GetCity()[i].GetY() << " ";
		output_old_square << path_square[best].GetCity()[i].GetAllele() << " " << path_square[best].GetDistance() << endl;
	}
	output_old_square << path_square[best].GetCity()[0].GetX() << " " << path_square[best].GetCity()[0].GetY() << " ";
	output_old_square << path_square[best].GetCity()[0].GetAllele() << " " << path_square[best].GetDistance() << endl;
	
	NewGeneration newgen_square(path_square);
	
	for(int i = 0; i < gen; i++) {
		newgen_square.Evolution(N_pop);
		newgen_square.Mutation();
	
		output_L_square << i << " " << newgen_square.GetPath()[0].GetDistance() << endl;
		for(int j = 0; j < N_pop / 2; j++) {
			dist += newgen_square.GetPath()[j].GetDistance();
		}
		output_meanL_square << i << " " << dist * 1. / (1.*N_pop/2.) << endl;
		dist = 0;
	}
	
	cout << "Final best path and distance on square " << endl;
	newgen_square.Fitness();
	newgen_square.GetPath()[best].PrintPath();
	newgen_square.GetPath()[best].PrintDistance();
	for(int i = 0; i < length_path; i++) {
		output_best_square << newgen_square.GetPath()[0].GetCity()[i].GetX() << " " << newgen_square.GetPath()[0].GetCity()[i].GetY() << " ";
		output_best_square << newgen_square.GetPath()[0].GetCity()[i].GetAllele() << " " << newgen_square.GetPath()[best].GetDistance() << " ";
		output_best_square << old_dist_square << endl;
	}
	output_best_square << newgen_square.GetPath()[0].GetCity()[0].GetX() << " " << newgen_square.GetPath()[0].GetCity()[0].GetY() << " ";	//Condizione di ritorno alla città iniziale
	output_best_square << newgen_square.GetPath()[0].GetCity()[0].GetAllele() << " " << newgen_square.GetPath()[best].GetDistance() << " ";
	output_best_square << old_dist_square << endl;
	
	output_L_square.close();
	output_meanL_square.close();
	output_best_square.close();
	output_old_square.close();
	
	return 0;
}
