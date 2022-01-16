CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio9.exe : Esercizio9.o GenAlg.o
	$(CC) GenAlg.o Esercizio9.o -o Esercizio9.exe
Esercizio9.o : Esercizio9.cpp GenAlg.h 
	$(CC) -c Esercizio9.cpp -o Esercizio9.o $(CFLAGS)
GenAlg.o : GenAlg.cpp GenAlg.h
	$(CC) -c GenAlg.cpp -o GenAlg.o $(CFLAGS)
clean :
	rm *.o Esercizio9.exe
