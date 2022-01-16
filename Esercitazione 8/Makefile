CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio8.exe : Esercizio8.o random.o
	$(CC) random.o Esercizio8.o -o Esercizio8.exe
Esercizio8.o : Esercizio8.cpp
	$(CC) -c Esercizio8.cpp -o Esercizio8.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio8.exe seed.out
