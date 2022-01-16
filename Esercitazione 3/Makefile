CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio3.exe : Esercizio3.o random.o
	$(CC) random.o Esercizio3.o -o Esercizio3.exe
Esercizio3.o : Esercizio3.cpp
	$(CC) -c Esercizio3.cpp -o Esercizio3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio3.exe seed.out
