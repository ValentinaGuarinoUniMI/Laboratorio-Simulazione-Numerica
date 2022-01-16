MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)
CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio10.exe : Esercizio10.o GenAlg.o random.o
	$(CC) $(MPI_COMPILE_FLAGS) random.o GenAlg.o Esercizio10.o  $(MPI_LINK_FLAGS) -o Esercizio10.exe
Esercizio10.o : Esercizio10.cpp GenAlg.h 
	$(CC) -c Esercizio10.cpp -o Esercizio10.o $(CFLAGS)
GenAlg.o : GenAlg.cpp GenAlg.h
	$(CC) -c GenAlg.cpp -o GenAlg.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio10.exe seed.out
