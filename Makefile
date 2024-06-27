OPT += -DEXTERNAL

CXX = icpx -Ofast -ffast-math -march=native -mavx2 -funroll-loops -std=c++17 -Wfatal-errors -fopenmp#-O3 -fp-model precise -march=core-avx2 -fopenmp -mavx -mavx2 -mfma #-ltcmalloc -mavx512fz
INCL += -I/usr/local/include -DH5_USE_16_API -I./ 
INCL += -L/usr/local/lib -lhdf5 -lz #-ltbb -D_GLIBCXX_USE_TBB_PAR_BACKEND=0


all:
	@make falcon

falcon: main.cxx
	$(CXX) $? -o $@ $(OPT) $(INCL)

clean:
	$(RM) ./*.o ./falcon ./fmm ./*~
