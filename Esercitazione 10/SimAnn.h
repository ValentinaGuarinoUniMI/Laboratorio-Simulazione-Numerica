#ifndef _GenAlg__h__
#define _GenAlg__h__

#include<stdio.h>
#include<cmath>
#include<vector>

using namespace std;
/************************************************************************************************************************************************************
	Classe Cities: classe per la formazione dei geni
************************************************************************************************************************************************************/
class Cities {								//Geni: ogni città rappresenta un gene che ha alleli diversi
private:
	int allele;
	double x, y;

public:
	Cities();					//Default constructor
	Cities(double, double, int);			//Costruttore con parametri
	Cities(const Cities&);				//Copy constructor
	Cities& operator=(const Cities&);		//Copy assignment
	bool operator==(const Cities&) const;
	bool operator!=(const Cities&) const;
	//~Cities(); 					//Trivial destructor
	void SetCoordinate(double, double);
	void SetAllele(int);

//Metodi per accedere alle variabili private

	double GetX() const;
	double GetY() const;
	int GetAllele() const;
	void PrintCity() const;
	
};

/************************************************************************************************************************************************************
	Classe Path: classe per la generazione del cromosoma
************************************************************************************************************************************************************/

class Path {		//Cromosoma: il percorso compiuto è dato dalla sequenza di più città

private:
	int length;
	double distance;
	vector <Cities> city;

public: 
				Path();									//Default constructor
	Path(vector<vector<double>>, int);		//Costruttore con parametri: prende in input la matrice formata da coordinate ed allele e la lunghezza
	Path(vector<Cities>);
	Path(const Path&);						//Copy constructor
	Path& operator=(const Path&);			//Copy assignment
	bool operator==(const Path&) const;
	bool operator!=(const Path&) const;
	//~Path(); 					//Trivial destructor
	void SetDistance();
	void SetLength(int);
	void PrintPath() const;
	void PrintDistance() const;
	double GetDistance() const;
	double GetLength() const;
	vector<Cities> GetCity()const;
	
	//Mutazioni
	void Pair_Permutation();
	void Contiguous_Permutation();
	void Inversion_Permutation();
	void Shift_Mutation();
	void Mutation();
	

};


void SetPosition_Circumference(vector<vector<double>>&, int);
void SetPosition_Square(vector<vector<double>>&, int);
void RandomPath(vector<vector<double>>&, int);
void Fitness(vector <Path>&);

#endif
