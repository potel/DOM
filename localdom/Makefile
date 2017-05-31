OBJECTS = scat.o potPara.o merson.o legendre.o waves.o minimize1D.o minimizeND.o fit.o imaginaryForm.o\
whit.o coul.o compound.o level.o channel.o expInt.o reaction.o volume.o surface.o hartreeFock.o spinOrbit.o\
 asy.o gaussInteg.o  gauss16.o sphericalB.o fitOM.o disperse.o twoFermi.o surfaceTF.o surVolume.o  levelD.o\
  compoundException.o standard.o surfaceSTD.o potplus.o calli.o non_local_wavefunction_3.o 
DEPS = test.h structs.h definiciones.h tremendo.h
HOLEDIR=./dom/ca48/gregory

CC=g++
LIBDIR=/usr/local/lib
CFLAGS=-I.:/usr/local/include -L$(LIBDIR) -lgsl -lgslcblas
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



CC=g++

###################
# Fortran Compiler# 
# and flag        #
###################
FC=gfortran 
FFLAGS= -c -O3 -fopenmp -fbounds-check -ffree-line-length-none
FLINK=gfortran

DIR := $(HOLEDIR)/build
DIR_RH7 := $(HOLEDIR)/build/rh7
DIR_COMMON_RH7 := ./dom/common/rh7

################
#  C++ files   #
################
SRC := spectralFunction.cpp partial.cpp fit.cpp reaction.cpp volume.cpp hartreeFock.cpp numerov.cpp pot.cpp scatterRspace.cpp boundRspace.cpp read_parameters.cpp compound.cpp plot.cpp

OBJ_RH7 := $(SRC:%.cpp=$(DIR_RH7)/%.o)

SRC_COMMON_RH7 := $(DIR_COMMON_RH7)/FourierBessel.cpp $(DIR_COMMON_RH7)/dom_nonlocal_rspace.cpp $(DIR_COMMON_RH7)/RedSelfEnergy.cpp $(DIR_COMMON_RH7)/kMesh.cpp $(DIR_COMMON_RH7)/clebsch_gordan.cpp $(DIR_COMMON_RH7)/potPara.cpp $(DIR_COMMON_RH7)/merson.cpp $(DIR_COMMON_RH7)/waves.cpp $(DIR_COMMON_RH7)/minimize1D.cpp $(DIR_COMMON_RH7)/minimizeND.cpp $(DIR_COMMON_RH7)/imaginaryForm.cpp $(DIR_COMMON_RH7)/whit.cpp $(DIR_COMMON_RH7)/coul.cpp $(DIR_COMMON_RH7)/level.cpp $(DIR_COMMON_RH7)/channel.cpp $(DIR_COMMON_RH7)/expInt.cpp $(DIR_COMMON_RH7)/spinOrbit.cpp $(DIR_COMMON_RH7)/asy.cpp $(DIR_COMMON_RH7)/gaussInteg.cpp $(DIR_COMMON_RH7)/gauss16.cpp $(DIR_COMMON_RH7)/disperse.cpp $(DIR_COMMON_RH7)/twoFermi.cpp $(DIR_COMMON_RH7)/surfaceTF.cpp $(DIR_COMMON_RH7)/surVolume.cpp $(DIR_COMMON_RH7)/levelD.cpp $(DIR_COMMON_RH7)/compoundException.cpp $(DIR_COMMON_RH7)/standard.cpp $(DIR_COMMON_RH7)/surfaceSTD.cpp $(DIR_COMMON_RH7)/numerovR.cpp $(DIR_COMMON_RH7)/io.cpp $(DIR_COMMON_RH7)/eigen.cpp $(DIR_COMMON_RH7)/density.cpp $(DIR_COMMON_RH7)/meshes.cpp $(DIR_COMMON_RH7)/numerical.cpp $(DIR_COMMON_RH7)/surfaceGeneral.cpp $(DIR_COMMON_RH7)/read_data.cpp $(DIR_COMMON_RH7)/sphericalB.cpp $(DIR_COMMON_RH7)/legendre.cpp

OBJ_COMMON_RH7 := $(SRC_COMMON_RH7:%.cpp=%.o)

#################
# Fortran files #
#################
SRC_FOR := subeep.f eep_dwia.f
OBJ_FOR := $(SRC_FOR:%.f=$(DIR)/%.o)
########################################

#Nothing above this block should be changed
##################################################################################################################
#These includes and libraries are machine-dependent, make sure they correspond to yours if you want this to work##
##################################################################################################################
# Includes
LOCAL_INCLUDE = -I./
COMMON_INCLUDE_RH7 = -I$(DIR_COMMON_RH7)
BNB_CFLAGS = -I/usr/local/numeric_bindings-v1
BOOST = -I/usr/local/boost_1_63_0
EIGEN_INCLUDE = -I/usr/local/eigen3
GSL_CFLAGS = $(shell gsl-config --cflags)

# Libraries
ATLAS_LIBS = -L/usr/local/Atlas/lib 
ALL_LIBS=-llapack -latlas  -lcblas    -llapack_atlas
GSL_LIBS = $(shell gsl-config --libs)
#################################################################################################################
#################################################################################################################

#This is the part that links the code to the libraries and includes all of the header files
CFLAGS= -c -O3  -Wfatal-errors  $(BNB_CFLAGS) $(EIGEN_INCLUDE) $(LOCAL_INCLUDE) $(GSL_CFLAGS) $(BOOST) -fopenmp -fpermissive
#LINKOPTION =   $(GSL_LIBS) $(ATLAS_LIBS) $(ALL_LIBS) 
LINKOPTION =   -llapack -lcblas  $(GSL_LIBS)

#################################################################################################################
#Nothing below this line should be changed

dphole: $(DIR_RH7)/myholegreen.o $(OBJ_RH7) $(OBJ_FOR) $(OBJ_COMMON_RH7) 
	$(CC)   -o dphole irreducible.o $(OBJECTS) $(DIR_RH7)/myholegreen.o  $(OBJ_COMMON_RH7) $(OBJ_FOR) $(OBJ_RH7) $(LINKOPTION) -l$(FLINK)  -fopenmp -I$(HOLEDIR)

$(DIR_RH7)/%.o: %.cpp | $(DIR_RH7)
	$(CC) $(CFLAGS) $(COMMON_INCLUDE_RH7) $< -o $@ 

$(DIR_COMMON_RH7)/%.o: $(DIR_COMMON_RH7)/%.cpp
	$(CC) $(CFLAGS) $(COMMON_INCLUDE_RH7) $< -o $@ 

$(DIR)/%.o : %.f | $(DIR)
	$(FC) $(FFLAGS) $< -o $@

$(DIR_RH7)/chisq_rh7.o: chisq.cpp | $(DIR_RH7)
	$(CC) $(CFLAGS) $(COMMON_INCLUDE_RH7) $< -o $@ 

$(DIR_RH7)/myholegreen.o: $(HOLEDIR)/myholegreen.cpp | $(DIR_RH7)
	$(CC) $(CFLAGS) $(COMMON_INCLUDE_RH7)  $< -o $@ -I$(HOLEDIR)
cleanhole:
	rm $(DIR)/*.o
cleanrh7:
	rm $(DIR_RH7)/*.o $(DIR)/*.o $(DIR_COMMON_RH7)/*.o *_rh7
cleanall:
	rm $(DIR)/*.o $(DIR_RH7)/*.o $(DIR_COMMON_RH7)/*.o *_rh6 *_rh7 holegreen

