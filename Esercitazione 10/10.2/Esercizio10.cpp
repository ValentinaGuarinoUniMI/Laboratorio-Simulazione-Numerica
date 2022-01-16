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
#include"GenAlg.h"
#include "mpi.h"
#include "random.h"
using namespace std;

int main(int argc, char* argv[]) {
	#define N 2
	MPI_Init(&argc, &argv);		//Inizializzazione della libreria di Message Passing

	int N_pop = 100;	//Grandezza della popolazione
	int length_path = 32;	//Numero di città
	int gen = 800;		//Numero di generazioni
	int size, rank;
	int tagA = 1;
	int tagB = 2;
	int tagC = 3;
	int tagD = 4;
	int best = 0;
	
	vector<vector<double>> position_square;
	vector <Path> path_square;
	int* alls_0 = new int[length_path];
	int* alls_1 = new int[length_path];
	int* alls_2 = new int[length_path];
	int* alls_3 = new int[length_path];
	
	int r_1 = 0, r_2 = 0, r_3 = 0, r_4 = 0;
	Random randgen;
	int seed[4]; 
	int p1, p2; 
	ifstream Primes("Primes"); 
	if (Primes.is_open()){   
		Primes >> p1 >> p2 ; 
	} else cerr << "PROBLEM: Unable to open Primes" << endl; 
	Primes.close(); 
	ifstream input("seed.in"); 
	string property; 
	if (input.is_open()){
		while ( !input.eof() ){ 
			input >> property; 
			if( property == "RANDOMSEED" ){ 
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3]; 
	 			randgen.SetRandom(seed,p1,p2); 
	 		} 
		} 
		input.close(); 
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat1, stat2, stat3, stat4;
	double tstart = MPI_Wtime();
//Ricerca percorso migliore per città poste all'interno di un quadrato

	//SetPosition_Square(position_square, length_path);
	SetPosition_FromFile_Square(position_square, length_path);	//Valuto l'andamento al variare delle generazioni del percorso 										  migliore iniziale di una data simulazione
	Path percorso(position_square, length_path);
	for (int i = 0; i < N_pop; i++) {
		if(i == 0) {
			path_square.push_back(Path(position_square, length_path));
		}
		else{
			RandomPath(position_square, length_path);	
			path_square.push_back(Path(position_square, length_path));
		}		
	}

	Fitness(path_square);

	
	if(rank == 0) {
		cout << "Initial best path and distance on square " << endl;
		path_square[best].PrintPath();
		path_square[best].PrintDistance();
	}
	
	NewGeneration newgen_square(path_square);
	
	for(int i = 0; i < gen; i++) {
		
		newgen_square.Evolution(N_pop);
		newgen_square.Mutation();
		if(i % 50 == 0) {			//Ogni 50 generazioni avviene lo scambio tra i nodi scelti in modo random
			r_1 = randgen.Rannyu_INT(0,3);
			r_2 = randgen.Rannyu_INT(0,3);
			r_3 = randgen.Rannyu_INT(0,3);
			r_4 = randgen.Rannyu_INT(0,3);
			
			while(r_1 == r_2) {
				r_2 = randgen.Rannyu_INT(0,3);	
			}
			while(r_2 == r_3 || r_3 == r_1) {
				r_3 = randgen.Rannyu_INT(0,3);	
			}
			while(r_3 == r_4 || r_4 == r_1 || r_4 == r_2) {
				r_4 = randgen.Rannyu_INT(0,3);	
			}
			//Inizializzazione dei vettori
			if(rank==r_1) {
				for(int j = 0; j < length_path; j++) {
					alls_0[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_1[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_2[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_3[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
				}
			}
			if(rank==r_2){
				for(int j = 0; j < length_path; j++) {
					alls_0[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_1[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_2[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_3[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
				}
			}
			if(rank==r_3){
				for(int j = 0; j < length_path; j++) {
					alls_0[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_1[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_2[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_3[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
				}
			}
			if(rank==r_4){
				for(int j = 0; j < length_path; j++) {
					alls_0[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_1[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_2[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
					alls_3[j] = newgen_square.GetPath()[best].GetCity()[j].GetAllele();
				}
			}
			//Scambio degli alleli tra coppie di nodi random	
			if(rank == r_1) {
				MPI_Send(&alls_0[0], length_path, MPI_INT, r_2, tagA, MPI_COMM_WORLD);
				MPI_Recv(&alls_1[0], length_path, MPI_INT, r_2, tagB, MPI_COMM_WORLD, &stat2);
			/*	cout << "Percorso scambiato: " << endl;
				for(int i = 0; i < length_path; i++) {
					cout << alls_0[i] << endl;
				}
				cout << "con " << endl;
				for(int i = 0; i < length_path; i++) {
					cout << alls_1[i] << endl;
				} */
				
				
			}
			
			if(rank == r_2) {
				MPI_Recv(&alls_0[0], length_path, MPI_INT, r_1, tagA, MPI_COMM_WORLD, &stat1);
				MPI_Send(&alls_1[0], length_path, MPI_INT, r_1, tagB, MPI_COMM_WORLD);
			}
			
			if(rank == r_3) {
				MPI_Send(&alls_2[0], length_path, MPI_INT, r_4, tagC, MPI_COMM_WORLD);
				MPI_Recv(&alls_3[0], length_path, MPI_INT, r_4, tagD, MPI_COMM_WORLD, &stat4);
			}
			if(rank == r_4) {
				MPI_Recv(&alls_2[0], length_path, MPI_INT, r_3, tagC, MPI_COMM_WORLD, &stat3);
				MPI_Send(&alls_3[0], length_path, MPI_INT, r_3, tagD, MPI_COMM_WORLD);
			}
			
			
			//Recupero delle coordinate per gli alleli scambiati
			if(rank == r_1) {
				Path new_path_0(percorso, alls_1, length_path);
				//cout << "Percorso scambiato " << endl;
				//new_path_0.PrintPath();
				//cout << "Miglior percorso vecchio " << endl;
				//newgen_square.GetPath()[best].PrintPath();
				newgen_square.SetPath(new_path_0, 0);
				//cout << "Miglior percorso nuovo per " << r_1 << endl;
				//newgen_square.GetPath()[best].PrintPath();
				
				
			}
			if(rank == r_2) {
				Path new_path_1(percorso, alls_0, length_path);
				newgen_square.SetPath(new_path_1, 0);
				//cout << "Miglior percorso nuovo per " << r_2 << endl;
				//newgen_square.GetPath()[best].PrintPath();
			}
			if(rank == r_3) {
				Path new_path_2(percorso, alls_3, length_path);
				newgen_square.SetPath(new_path_2, 0);
				//cout << "Miglior percorso nuovo per " << r_3 << endl;
				//newgen_square.GetPath()[best].PrintPath();
			}
			if(rank == r_4) {
				Path new_path_3(percorso, alls_2, length_path);
				newgen_square.SetPath(new_path_3, 0);
				//cout << "Miglior percorso nuovo per " << r_4 << endl;
				//newgen_square.GetPath()[best].PrintPath();
			}
						
		}
		
		Print(newgen_square, i, length_path, N_pop, rank);	//WARNING: file L e meanL salvati in append
		
	}
	
	cout << "Final best path and distance on square for node " << rank << endl;
	newgen_square.Fitness();
	newgen_square.GetPath()[best].PrintPath();
	cout << "Distance: ";
	newgen_square.GetPath()[best].PrintDistance();
	cout << endl;
	double tend = MPI_Wtime();
	double dt = tend - tstart;
    	cout << "Time measured: " << dt << " s" << endl;		//Calcolo del tempo di esecuzione
    	MPI_Finalize();
	return 0;
}
