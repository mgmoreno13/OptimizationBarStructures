# Makefile for compiling Meca
# Programa para an‡lisis matricula de estructuras

# Copyright (C) 2011, by Salvador Botello Rionda
# botello@cimat.mx
# Centro de Investigaci'on en Matem'aticas A.C.
CC=g++
RM=rm -f *~ *~*
CFLAGS= -O2 -fopenmp
#CFLAGS=-fopenmp -g
OBJS:= datos.o efi_CR.o efi_RC.o efi_RF_c.o efi_RF_o.o eficiencias.o estatico.o fuerzas.o  memoria.o pospro.o posprogid.o prelim.o principal.o raros.o rigidez.o solver.o tensiones.o  Optimizador.o valorespropios.o 

MAIN=meca

all:$(MAIN)
$(MAIN):$(OBJS) main.o
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) main.o
datos.o:  datos.cpp
	$(CC) $(CFLAGS) -c datos.cpp
efi_CR.o : efi_CR.cpp 
	$(CC) $(CFLAGS) -c efi_CR.cpp 
efi_RC.o : efi_RC.cpp
	$(CC) $(CFLAGS) -c efi_RC.cpp 
efi_RF_c.o : efi_RF_c.cpp
	$(CC) $(CFLAGS) -c efi_RF_c.cpp 
efi_RF_o.o : efi_RF_o.cpp
	$(CC) $(CFLAGS) -c efi_RF_o.cpp 
eficiencias.o : eficiencias.cpp
	$(CC) $(CFLAGS) -c eficiencias.cpp 
estatico.o : estatico.cpp
	$(CC) $(CFLAGS) -c estatico.cpp 
fuerzas.o : fuerzas.cpp
	$(CC) $(CFLAGS) -c fuerzas.cpp 
main.o : main.cpp
	$(CC) $(CFLAGS) -c main.cpp 
memoria.o : memoria.cpp
	$(CC) $(CFLAGS) -c memoria.cpp 
pospro.o : pospro.cpp
	$(CC) $(CFLAGS) -c pospro.cpp 
posprogid.o : posprogid.cpp
	$(CC) $(CFLAGS) -c posprogid.cpp 
prelim.o : prelim.cpp
	$(CC) $(CFLAGS) -c prelim.cpp 
principal.o : principal.cpp
	$(CC) $(CFLAGS) -c principal.cpp 
raros.o : raros.cpp
	$(CC) $(CFLAGS) -c raros.cpp 
rigidez.o : rigidez.cpp
	$(CC) $(CFLAGS) -c rigidez.cpp 
solver.o :solver.cpp
	$(CC) $(CFLAGS) -c solver.cpp 
tensiones.o : tensiones.cpp
	$(CC) $(CFLAGS) -c tensiones.cpp 
valorespropios.o : valorespropios.cpp
	$(CC) $(CFLAGS) -c valorespropios.cpp 
optimizador.o : Optimizador.cpp
	$(CC) $(CFLAGS) -c Optimizador.cpp

clean:
	$(RM) $(OBJS)  $(OBJMPI) meca.o 
cleanall:
	$(RM) $(OBJS) $(OBJMPI)  meca.o meca


