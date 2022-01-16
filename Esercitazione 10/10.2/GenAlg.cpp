#define _USE_MATH_DEFINES
#include<iostream>
#include<stdio.h>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<vector>
#include<algorithm>
#include<fstream>
#include<random>

#include "GenAlg.h"

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


double Cities::GetAllele() const {

	return allele;
}

/*void Cities::ShiftAllele(int* alls) {
	
	allele = &alls;

}*/

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


Path::Path(Path& percorso, int* allele, int len) {			
	
	length = len;
	city.resize(len);
	for(int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {			
			if(allele[i] == percorso.GetCity()[j].GetAllele()) {
				city[i].SetCoordinate(percorso.GetCity()[j].GetX(), percorso.GetCity()[j].GetY());
				city[i].SetAllele(percorso.GetCity()[j].GetAllele());
				break;
			}
			else {
				continue;
			}		
		}
	}
	
	SetDistance();
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


/************************************************************************************************************************************************************
	Classe New Generation
************************************************************************************************************************************************************/

NewGeneration::NewGeneration() {
//Deafault constructor
	new_length_path = 0;
	distance = 0;
}


NewGeneration::NewGeneration(vector<Path>& path) {
//Inizializzazione costruttore con parametri
	new_path = path;
	new_length_path = path.size();
}

NewGeneration::NewGeneration(const NewGeneration& original_class) {	//Copy constructor
	new_length_path = original_class.GetLengthPop();
	new_path = original_class.GetPath();
}

NewGeneration& NewGeneration::operator= (const NewGeneration& original_class) {	//Copy assignment
	new_length_path = original_class.GetLengthPop();
	new_path = original_class.GetPath();
	return *this;
}

void NewGeneration::SetLengthPop(int len) {
	new_length_path = len;
}

double NewGeneration::GetLengthPop() const {

	return new_length_path;
}


int NewGeneration::Selector(int length_pop) {	//Metodo per la selezione 
	double r, p = 0;
	int j = 0;
	
	
	random_device ran;
	default_random_engine gen(ran());
	uniform_real_distribution<> dist_r(0.0, 1.0);

	p = 2.5;
	for(int i = 0; i < length_pop; i++) {

		r = dist_r(gen);
		j = int(length_pop * pow(r, p));
		
	}
	
	return j;

}

vector<Path> NewGeneration::GetPath() const {
	return new_path;
}

void NewGeneration::SetPath(Path& percorso, int n) {
	new_path[n] = percorso;
}

void NewGeneration::Crossover(int mum, int dad){
	vector<Cities> parent_1 = new_path[mum].GetCity();
	vector<Cities> parent_2 = new_path[dad].GetCity();
	vector<Cities> son_1 = new_path[mum].GetCity();		//Inizialmente i due figli sono copie dei genitori
	vector<Cities> son_2 = new_path[dad].GetCity();
	
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_distr(1, new_path[mum].GetLength()-2);
	
	int cut = unif_distr(gen);
	//cout << "cut " << cut << endl;
	son_1.erase(son_1.begin() + cut, son_1.end());		//Si elimina un pezzo di percorso di un vettore figlio partendo da un 
	son_2.erase(son_2.begin() + cut, son_2.end());		//punto indicato dal numero casuale cut
	
	/*cout << "Primo genitore\n" << endl;
	for(int i = 0; i < parent_1.size(); i++){
		parent_1[i].PrintCity();
	}
	
	cout << "Secondo genitore\n" << endl;
	for(int i = 0; i < parent_2.size(); i++){
		parent_2[i].PrintCity();
	}
	*/
	
	//Crossover figlio 1
	for(int i = 1; i < new_path[dad].GetLength(); i++){
		if(find(son_1.begin(), son_1.end(), parent_2[i]) != son_1.end()){	//Viene attuato il crossover riempiendo la parte 
			continue;							//tagliata del vettore figlio 1 con gli elementi
		}									//contenuti nel vettore genitore 2 che non si 
		else {								 	//trovano già in figlio 1, mantenendo l'ordine
			son_1.push_back(parent_2[i]);					//in cui si trovano in genitore 2.
		}
	}
	
	//Crossover figlio 2
	for(int i = 1; i < new_path[mum].GetLength(); i++){
		if(find(son_2.begin(), son_2.end(), parent_1[i]) != son_2.end()){
			continue;
		}
		else {
			son_2.push_back(parent_1[i]);
		}
	
	}
	
	/*cout << "Primo figlio\n" << endl;
	for(int i = 0; i < son_1.size(); i++){
		son_1[i].PrintCity();
	} 
	
	cout << "Secondo figlio\n" << endl;
	for(int i = 0; i < son_2.size(); i++){
		son_2[i].PrintCity();
	}*/
	
	Path my_first_son(son_1);
	Path my_second_son(son_2);
	
	new_path.push_back(my_first_son);
	new_path.push_back(my_second_son);
}

void NewGeneration::Pair_Permutation(int index) {
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_swap(1,new_path[0].GetLength()-2);
	
	int swap_index = unif_swap(gen);
	vector<Cities> city_swapped = new_path[index].GetCity();
	
	/*cout << "Before permutation\n" << endl;
	
	for(int i = 0; i < city_swapped.size(); i++){
		city_swapped[i].PrintCity();
	}*/
	
	if(swap_index != new_path[0].GetLength()-2) {
		swap(city_swapped[swap_index], city_swapped[swap_index + 1]);
	}
	else {
		swap(city_swapped[swap_index], city_swapped[swap_index - 1]);
	}
	
	/*cout << "After permutation\n" << endl;
	
	for(int i = 0; i < city_swapped.size(); i++){
		city_swapped[i].PrintCity();
	}*/
	Path permutation(city_swapped);
	permutation.SetDistance();
	new_path[index] = permutation;

}


void NewGeneration::Inversion_Permutation(int index) {
	int first_ind = 0, second_ind = 0;
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_distr(1, new_path[0].GetLength() - 1);
	
	first_ind = unif_distr(gen);
	second_ind = unif_distr(gen);
	while(first_ind >= second_ind) {
		first_ind = unif_distr(gen);
		second_ind = unif_distr(gen);
	}
	//cout << "First index " << first_ind << " second index " << second_ind << endl;
	
	vector<Cities> inverted_path = new_path[index].GetCity();
	/*cout << "Before inversion\n" << endl;
	
	for(int i = 0; i < inverted_path.size(); i++){
		inverted_path[i].PrintCity();
	}*/
	
	
	reverse(inverted_path.begin()+first_ind, inverted_path.begin()+ second_ind + 1);
	
	/*cout << "After inversion\n" << endl;
	
	for(int i = 0; i < inverted_path.size(); i++){
		inverted_path[i].PrintCity();
	}*/
	
	Path inversion(inverted_path);
	inversion.SetDistance();
	new_path[index] = inversion;


}

void NewGeneration::Contiguous_Permutation(int index) {
	int m = 0;
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_distr(1, new_path[0].GetLength() / 2 - 1);
	
	m = unif_distr(gen);
	vector<Cities> cont_path = new_path[index].GetCity();
	rotate(cont_path.begin() + 1, cont_path.begin() + m, cont_path.begin() + 2*m);
	
	Path contig_perm(cont_path);
	contig_perm.SetDistance();
	new_path[index] = contig_perm;


}

void NewGeneration::Shift_Mutation(int index) {
	int m = 0;
	random_device rand;
	default_random_engine gen(rand());
	uniform_int_distribution<> unif_distr(1, new_path[0].GetLength() / 2 - 1);
	
	m = unif_distr(gen);
	vector<Cities> shifted_path = new_path[index].GetCity();
	rotate(shifted_path.begin() + 1, shifted_path.begin() + m, shifted_path.end());
	
	Path shifted(shifted_path);
	shifted.SetDistance();
	new_path[index] = shifted;


}

void NewGeneration::Evolution(int len) {
	int mum, dad = 0;
	double r = 0;
	int old_size = new_path.size();
	
	random_device ran;
	default_random_engine gen(ran());
	uniform_real_distribution<> dist_r(0.0, 1.0);
	
	for(int i = 0; i < len / 2.; i++) {
		r = dist_r(gen);	
		mum = Selector(len);
		dad = Selector(len);
		while (dad == mum) {
			dad = Selector(len);
		}
		//cout << mum << " " << dad << endl;
		if(r < 0.7) {
			
			//cout<<"Chiamo crossover tra "<<mum<<" "<<dad<<endl;
			Crossover(mum, dad);
		}
		else {
			
			//cout << "Niente crossover per "<<mum<< " " <<dad<<endl;
			new_path.push_back(new_path[mum]);
			new_path.push_back(new_path[dad]);
		}
		
	}
	
	new_path.erase(new_path.begin(), new_path.begin()+old_size);

	/*for (auto i : path) {
		i.PrintPath();
		cout << "\n" << endl;
	      //i.PrintDistance();
	}*/
	
}


void NewGeneration::Mutation() {
	random_device ran;
	default_random_engine gen(ran());
	uniform_real_distribution<> unif_random(0.0, 1.0);
	uniform_int_distribution<> unif_index(0, GetLengthPop() - 1);
	
	double first_coin, second_coin, third_coin, fourth_coin = 0;
	int index = 0;
	first_coin = unif_random(gen);
	second_coin = unif_random(gen);
	third_coin = unif_random(gen);
	fourth_coin = unif_random(gen);
	index = unif_index(gen);
	
	for(int i = 0; i < new_length_path; i++) {
	
		//Ogni mutazione ha probabilità del 10% di avvenire
		if(first_coin < 0.1) {
			//cout << "Pair permutation on " << index << endl;
			Pair_Permutation(index);
		}
		index = unif_index(gen);
		if(second_coin < 0.1) {
			//cout << "Inversion permutation on "<< index << endl;
			Inversion_Permutation(index);
		}
		index = unif_index(gen);
		if(third_coin < 0.1) {
			Contiguous_Permutation(index);
		}
		index = unif_index(gen);
		if(fourth_coin < 0.1) {
			Shift_Mutation(index);
		}	
	}
	
	Fitness();
}


void NewGeneration::Fitness()
{
	sort(new_path.begin(), new_path.end(), [](const Path& lhs, const Path& rhs) {	//Lambda function per il sorting
		return lhs.GetDistance() < rhs.GetDistance();
		});
}


/************************************************************************************************************************************************************
	Funzioni
************************************************************************************************************************************************************/
void SetPosition_Circumference(vector<vector<double>>& pos, int len_path) {
	double theta;
	double r = 10.0;
	random_device ran;
	default_random_engine gen(ran());
	uniform_real_distribution<> dist(0., 2. * M_PI);
	
	for (int i = 0; i < len_path; i++) {
		theta = dist(gen);
		pos.push_back({ r * cos(theta), r * sin(theta), (i + 1) * 1. });
	}

}

void SetPosition_Square(vector<vector<double>>& pos, int len_path) {
	double x = 0, y = 0;
	
	random_device ran;
	default_random_engine gen(ran());
	uniform_real_distribution<> dist(-1, 1);
	
	for (int i = 0; i < len_path; i++) {
		x = dist(gen);
		y = dist(gen);
		pos.push_back({ x, y, (i + 1) * 1. });
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


void SetPosition_FromFile(vector<vector<double>>& pos, int len_path) {	//Funzione modificata per prendere in input come percorso iniziale quello corrispondente al miglior percorso iniziale di una data simulazione e valutarne l'andamento al variare delle generazioni
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

void SetPosition_FromFile_Square(vector<vector<double>>& pos, int len_path) {	//Funzione modificata per prendere in input come percorso iniziale quello corrispondente al miglior percorso iniziale di una data simulazione e valutarne l'andamento al variare delle generazioni
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

void Print(NewGeneration& newgen, int gen, int len_path, int N_pop, int rank) {

ofstream out_best;
ofstream out_L;
ofstream out_meanL;
double dist = 0;
int best = 0;

if(rank == 0) {
	out_best.open("bestpath_node0.dat");
	out_L.open("L_node0.dat",ios::app);
	out_meanL.open("meanL_node0.dat", ios::app);
}

else if(rank == 1) {
	out_best.open("bestpath_node1.dat");
	out_L.open("L_node1.dat", ios::app);
	out_meanL.open("meanL_node1.dat", ios::app);
}
  
else if(rank == 2) {
	out_best.open("bestpath_node2.dat");
	out_L.open("L_node2.dat", ios::app);
	out_meanL.open("meanL_node2.dat", ios::app);
}	

else {
	out_best.open("bestpath_node3.dat");
	out_L.open("L_node3.dat", ios::app);
	out_meanL.open("meanL_node3.dat", ios::app);
}

out_L << gen << " " << newgen.GetPath()[0].GetDistance() << endl;

for(int j = 0; j < N_pop / 2; j++) {
		dist += newgen.GetPath()[j].GetDistance();
}
out_meanL << gen << " " << dist * 1. / (1.* N_pop/2.) << endl;
dist = 0;

for(int i = 0; i < len_path; i++) {
		out_best << newgen.GetPath()[0].GetCity()[i].GetX() << " " << newgen.GetPath()[0].GetCity()[i].GetY() << " ";
		out_best << newgen.GetPath()[0].GetCity()[i].GetAllele() << " " << newgen.GetPath()[best].GetDistance() << " " << endl;
	}
out_best << newgen.GetPath()[0].GetCity()[0].GetX() << " " << newgen.GetPath()[0].GetCity()[0].GetY() << " ";	//Condizione di ritorno alla città iniziale
out_best << newgen.GetPath()[0].GetCity()[0].GetAllele() << " " << newgen.GetPath()[best].GetDistance() << " "<< endl;
  	
}
