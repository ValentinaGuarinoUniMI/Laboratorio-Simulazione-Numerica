CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio5.exe : Esercizio5.o random.o
	$(CC) random.o Esercizio5.o -o Esercizio5.exe
Esercizio5.o : Esercizio5.cpp
	$(CC) -c Esercizio5.cpp -o Esercizio5.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio5.exe seed.out
