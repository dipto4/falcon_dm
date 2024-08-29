#enable this for additional debugging
#OPT += -DDEBUG


CXX = icpx # change it to your favorite compiler
CXXFLAGS = -std=c++17 -Ofast -ffast-math -march=native -mavx2 -funroll-loops -Wfatal-errors -fopenmp

FALCONDIR=./
H5DIR=/usr/local/include
BOOSTDIR=./

H5LIB=/usr/local/lib

INCL += -I$(H5DIR) -DH5_USE_16_API -I$(BOOSTDIR) -I$(FALCONDIR) 
INCL += -L$(H5LIB) -lhdf5 -lz #-ltbb -D_GLIBCXX_USE_TBB_PAR_BACKEND=0


all:
	@make falcon

falcon: main.cxx
	$(CXX) $(CXXFLAGS) $? -o $@ $(OPT) $(INCL)

clean:
	$(RM) ./*.o ./falcon 
