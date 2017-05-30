OBJECTS = scat.o potPara.o merson.o legendre.o waves.o minimize1D.o minimizeND.o fit.o imaginaryForm.o\
whit.o coul.o compound.o level.o channel.o expInt.o reaction.o volume.o surface.o hartreeFock.o spinOrbit.o\
 asy.o gaussInteg.o  gauss16.o sphericalB.o fitOM.o disperse.o twoFermi.o surfaceTF.o surVolume.o  levelD.o\
  compoundException.o standard.o surfaceSTD.o potplus.o calli.o  non_local_wavefunction_3.o
DEPS = test.h structs.h definiciones.h tremendo.h guass_wg.h
CC=g++
LIBDIR=/home/ansheng/gsl/lib
CFLAGS=-I.:/home/ansheng/gsl/include -L$(LIBDIR) -lgsl -lgslcblas -larmadillo
#CFLAGS=-L$(LIBRARY_PATH)   -lgsl -lgslcblas
#CFLAGS=-L$(LIBRARY_PATH) -L/opt/software/GCC/4.4.5/lib -I/opt/software/GCC/4.4.5/include  -lgsl -lgslcblas
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
#.f.o:
#	${FC} -c -fbounds-check subeep.f
dpdom: irreducible.o $(OBJECTS) 
	g++ -o dpdom  $(OBJECTS) $(CPPFLAGS) $(CFLAGS) 

clean:
	rm -f *.o irreducible

