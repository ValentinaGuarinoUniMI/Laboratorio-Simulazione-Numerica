CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio10.exe : Esercizio10.o SimAnn.o
	$(CC) SimAnn.o Esercizio10.o -o Esercizio10.exe
Esercizio10.o : Esercizio10.cpp SimAnn.h 
	$(CC) -c Esercizio10.cpp -o Esercizio10.o $(CFLAGS)
SimAnn.o : SimAnn.cpp SimAnn.h
	$(CC) -c SimAnn.cpp -o SimAnn.o $(CFLAGS)
clean :
	rm *.o Esercizio10.exe
