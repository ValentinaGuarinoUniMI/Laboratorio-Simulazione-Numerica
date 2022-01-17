#define _USE_MATH_DEFINES
#include<iostream>
#include<stdio.h>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<vector>
#include<algorithm>
#include<fstream>
#include <iterator>  // for istream_iterator
#include<random>
#include <string> 
#include <sstream>

#include "SimAnn.h"

using namespace std;

/************************************************************************************************************************************************************
	Classe Cities
************************************************************************************************************************************************************/

Cities::Cities() {
//Inizializzazione default constructor
	x = 0;
	y = 0;
	allele = 0;
}

Cities::Cities(double coord_x, double coord_y, int all) { 
//Inizializzazione costruttore con parametri
	x = coord_x;
	y = coord_y;
	allele = all;

}

Cities::Cities(const Cities& original_class) { //Copy constructor
	x = original_class.GetX();
	y = original_class.GetY();
	allele = original_class.GetAllele();
}

Cities& Cities::operator= (const Cities &original_class) { //Copy assignment
	x = original_class.GetX();
	y = original_class.GetY();
	allele = original_class.GetAllele();
	return *this;
}

bool Cities::operator==(const Cities& is_equal) const { //Overload operator ==
	if((x==is_equal.GetX()) && 
	   (y==is_equal.GetY()) &&
	   (allele==is_equal.GetAllele()))
	   return true;
	else
	   return false;
}

bool Cities::operator!=(const Cities& is_equal) const {
//Overload operator !=
	if((x!=is_equal.GetX()) && 
	   (y!=is_equal.GetY()) &&
	   (allele!=is_equal.GetAllele()))
	   return true;
	else
	   return false;
}


void Cities::SetCoordinate(double coord_x, double coord_y) {
	
	x = coord_x;
	y = coord_y;
}

void Cities::SetAllele(int all) {

	allele = all;

}

double Cities::GetX() const {

	return x;
}

double Cities::GetY() const {

	return y;
}


int Cities::GetAllele() const {

	return allele;
}

void Cities::PrintCity() const {

	cout << x << " " << y << " " << allele << endl;
}


/************************************************************************************************************************************************************
	Classe Path
************************************************************************************************************************************************************/

Path::Path() {						//Inizializzazione default constructor
	length = 0;
	distance = 0;
}


Path::Path(vector<vector<double>> pos, int len) {			//Inizializzazione costruttore con parametri

	length = len;
	city.resize(len);
	
	for (int i = 0; i < city.size(); i++) {			//i: righe della matrice(x, y e allele per le diverse città)
		int j = 0;									//j: colonne della matrice (0: x, 1: y, 2: allele)
		city[i].SetCoordinate(pos[i][j], pos[i][j + 1]);
		city[i].SetAllele(pos[i][j + 2]);
	}
	SetDistance();
}

Path::Path(vector<Cities> my_city) {
//Secondo costruttore con parametri
	city = my_city;
	length = my_city.size();
	double dist = 0;
	for (int i = 0; i < my_city.size(); i++) {
		int j = i + 1;
		double x_dist, y_dist;
		if (i == length - 1) {							//Condizione di ritorno alla prima città
			x_dist = my_city[i].GetX() - my_city[0].GetX();
			y_dist = my_city[i].GetY() - my_city[0].GetY();

			dist += sqrt(pow(x_dist,2) + pow(y_dist,2));		//Norma L
		}
		else {
			x_dist = my_city[i].GetX() - my_city[j].GetX();
			y_dist = my_city[i].GetY() - my_city[j].GetY();
	
			dist += sqrt(pow(x_dist, 2) + pow(y_dist, 2));	//Norma L
		}
	}
	distance = dist;

}

Path::Path(const Path& original_class) {	//Copy constructor
	length = original_class.GetLength();
	distance = original_class.GetDistance();
	city = original_class.GetCity();
}

Path& Path::operator= (const Path& original_class) {	//Copy assignment
	length = original_class.GetLength();
	distance = original_class.GetDistance();
	city = original_class.GetCity();
	return *this;
}

bool Path::operator==(const Path& is_equal) const{	//Overload operator ==
	if((length==is_equal.GetLength()) &&
	   (distance==is_equal.GetDistance()) && 
	   (city==is_equal.GetCity()))
	   return true;
 	else
 	   return false;	
}

bool Path::operator!=(const Path& is_equal) const{	//Overload operator !=
	if((length!=is_equal.GetLength()) &&
	   (distance!=is_equal.GetDistance()) && 
	   (city!=is_equal.GetCity()))
	   return true;
 	else
 	   return false;	
}

void Path::SetLength(int len) {
	length = len;
}

void Path::SetDistance() {							//N.B.: trovandosi le città su una circonferenza, alcuni percorsi
												    // potrebbero condividere la stessa distanza per simmetria.

	double dist = 0;
	for (int i = 0; i < length; i++) {
		int j = i + 1;
		double x_dist, y_dist;
		if (i == length - 1) {							//Condizione di ritorno alla prima città
			x_dist = city[i].GetX() - city[0].GetX();
			y_dist = city[i].GetY() - city[0].GetY();

			dist += sqrt(pow(x_dist,2) + pow(y_dist,2));		//Norma L
		}
		else {
			x_dist = city[i].GetX() - city[j].GetX();
			y_dist = city[i].GetY() - city[j].GetY();
	
			dist += sqrt(pow(x_dist, 2) + pow(y_dist, 2));	//Norma L
		}
	}
	distance = dist;
}

double Path::GetDistance() const {

	return distance;
}

double Path::GetLength() const {

	return length;
}

vector<Cities> Path::GetCity() const {
	return city;
}

void Path::PrintDistance() const {

	GetDistance();
	cout << distance << endl;
}

void Path::PrintPath() const {
	for (auto i : city) {
		i.PrintCity();
	}
}

void Path::Pair_Permutation() {
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_swap(1,length - 2);
	
	int swap_index = unif_swap(gen);
	vector<Cities> city_swapped = city;
	
	/*cout << "Before permutation\n" << endl;
	
	for(int i = 0; i < city_swapped.size(); i++){
		city_swapped[i].PrintCity();
	}*/
	
	if(swap_index != length-2) {
		swap(city_swapped[swap_index], city_swapped[swap_index + 1]);
	}
	else {
		swap(city_swapped[swap_index], city_swapped[swap_index - 1]);
	}
	
	/*cout << "After permutation\n" << endl;
	
	for(int i = 0; i < city_swapped.size(); i++){
		city_swapped[i].PrintCity();
	}*/
	
	city = city_swapped;
	SetDistance();

}


void Path::Inversion_Permutation() {
	int first_ind = 0, second_ind = 0;
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_distr(1, length - 1);
	
	first_ind = unif_distr(gen);
	second_ind = unif_distr(gen);
	while(first_ind >= second_ind) {
		first_ind = unif_distr(gen);
		second_ind = unif_distr(gen);
	}
	//cout << "First index " << first_ind << " second index " << second_ind << endl;
	
	vector<Cities> inverted_path = city;
	/*cout << "Before inversion\n" << endl;
	
	for(int i = 0; i < inverted_path.size(); i++){
		inverted_path[i].PrintCity();
	}*/
	
	
	reverse(inverted_path.begin()+first_ind, inverted_path.begin()+ second_ind + 1);
	
	/*cout << "After inversion\n" << endl;
	
	for(int i = 0; i < inverted_path.size(); i++){
		inverted_path[i].PrintCity();
	}*/
	
	city = inverted_path;
	SetDistance();

}

void Path::Contiguous_Permutation() {
	int m = 0;
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_distr(1, length / 2 - 1);
	
	m = unif_distr(gen);
	vector<Cities> cont_path = city;
	rotate(cont_path.begin() + 1, cont_path.begin() + m, cont_path.begin() + 2*m);
	
	city = cont_path;
	SetDistance();
}

void Path::Shift_Mutation() {
	int m = 0;
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_distr(1, length / 2 - 1);
	
	m = unif_distr(gen);
	vector<Cities> shifted_path = city;
	rotate(shifted_path.begin() + 1, shifted_path.begin() + m, shifted_path.end());
	
	city = shifted_path;
	SetDistance();

}

void Path::Mutation() {
	
	Pair_Permutation();
	Inversion_Permutation();
	Contiguous_Permutation();
	Shift_Mutation();
	

}
/************************************************************************************************************************************************************
	Funzioni
************************************************************************************************************************************************************/
void SetPosition_Circumference(vector<vector<double>>& pos, int len_path) {	//Funzione modificata per prendere in input come percorso iniziale quello corrispondente al miglior percorso della prima generazione dell'Esercizio 9, in modo da poter effettuare il confronto

	ifstream input;
	input.open("initial_path.dat");
	
	vector<double> positions;
	double a = 0, b = 0;
		
	while(b != len_path)
 	{
 		while(positions.size() != 3)
 		{
 			input >> a;
 			positions.push_back(a);
 		}
 	
 	
 		pos.push_back(positions);
 		positions.clear();
 		b++;
 	}
  	
}

void SetPosition_Square(vector<vector<double>>& pos, int len_path) { //Funzione modificata per prendere in input come percorso iniziale quello corrispondente al miglior percorso della prima generazione dell'Esercizio 9, in modo da poter effettuare il confronto
	ifstream input;
	input.open("initial_path_square.dat");
	
	vector<double> positions;
	double a = 0, b = 0;
		
	while(b != len_path)
 	{
 		while(positions.size() != 3)
 		{
 			input >> a;
 			positions.push_back(a);
 		}
 	
 	
 		pos.push_back(positions);
 		positions.clear();
 		b++;
 	}

}


void RandomPath(vector<vector<double>>& pos, int len_path) {

	random_shuffle(pos.begin() + 1, pos.end());	//begin + 1 perché la città di partenza è sempre la stessa
}												//N.B.: poiché la prima città è fissata, scegliendo 5 città iniziali 
												//ho solo 4!=24 possibili percorsi diversi. Prendendo 32 possibili percorsi ne avremo allora alcuni ripetuti.

void Fitness(vector<Path>& path)
{
	sort(path.begin(), path.end(), [](const Path& lhs, const Path& rhs) {	//Lambda function per il sorting
		return lhs.GetDistance() < rhs.GetDistance();
		});
}

