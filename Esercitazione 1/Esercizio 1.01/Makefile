CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio1.exe : Esercizio1.o random.o
	$(CC) random.o Esercizio1.o -o Esercizio1.exe
Esercizio1.o : Esercizio1.cpp
	$(CC) -c Esercizio1.cpp -o Esercizio1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio1.exe seed.out