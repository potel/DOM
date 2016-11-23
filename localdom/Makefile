OBJECTS = scat.o potPara.o merson.o legendre.o waves.o minimize1D.o minimizeND.o fit.o imaginaryForm.o\
whit.o coul.o compound.o level.o channel.o expInt.o reaction.o volume.o surface.o hartreeFock.o spinOrbit.o\
 asy.o gaussInteg.o  gauss16.o sphericalB.o fitOM.o disperse.o twoFermi.o surfaceTF.o surVolume.o  levelD.o\
  compoundException.o standard.o surfaceSTD.o potplus.o calli.o 
#OBJECTS =    potplus.o calli.o reaction.o 

#OBJECTS2 =  scat.o PotPara.o merson.o Legendre.o waves.o whit.o coul.o ExpInt.o
#ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
#DEPS = test.h structs.h definiciones.h tremendo.h
CC=g++
LIBDIR=/usr/local/lib
CFLAGSX=-I.:/usr/local/include -L$(LIBDIR) -lgsl -lgslcblas
CFLAGS=-L$(LIBRARY_PATH)   -lgsl -lgslcblas
#CFLAGS=-L$(LIBRARY_PATH) -L/opt/software/GCC/4.4.5/lib -I/opt/software/GCC/4.4.5/include  -lgsl -lgslcblas
#%.o: %.c $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS)
#.f.o:
#	${FC} -c -fbounds-check subeep.f
dpdom: irreducible.o $(OBJECTS) 
	g++ -o dpdom  $(OBJECTS) $(CPPFLAGS) $(CFLAGS) 
#dpdomx: irreducible.o $(OBJECTS) 
	$(CC) -o dpdom  $(OBJECTS)  $(CFLAGSX) 

#irreducible: irreducible.o $(OBJECTS) 
#	$(CC) -o irreducible irreducible.o $(OBJECTS) $(CPPFLAGS) $(CFLAGS) 
#irreducible: irreducible.o $(OBJECTS) 
#	$(CC) -o irreducible irreducible.o $(OBJECTS) $(CFLAGS) -lgsl -lgslcblas


#$(ALLOBJECTS): %.o : %.cpp 
#	$(CC) $(CFLAGS) $< -o $@


clean:
	rm -f *.o irreducible

