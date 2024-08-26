OPT += -DEXTERNAL

CXX = icpx -Ofast -ffast-math -march=native -mavx2 -funroll-loops -std=c++17 -Wfatal-errors -fopenmp#-O3 -fp-model precise -march=core-avx2 -fopenmp -mavx -mavx2 -mfma #-ltcmalloc -mavx512fz
#CXX = icpx -std=c++17 -funroll-loops -O3 -fp-model precise -march=broadwell -qopenmp -fma #-mavx2#-fomit-frame-pointer #-ltcmalloc -mavx512fz
#CXX = icpc -std=c++17 -funroll-loops -O3 -fp-model precise -march=native -fopenmp -fcx-limited-range -mavx -mavx2 -mfma #-ltcmalloc -mavx512f
#CXX = g++ -funroll-loops -Wfatal-errors -O3 -Wno-format -mavx -march=native -fopenmp -mavx -mavx2 -mfma
#INCL += -I/hildafs/home/diptajym/hildafs/dependencies/hdf5_clang_install/include -DH5_USE_16_API -I./
#INCL += -L/hildafs/home/diptajym/hildafs/dependencies/hdf5_clang_install/lib -lhdf5 -lz -ltbb -D_GLIBCXX_USE_TBB_PAR_BACKEND=0  

INCL += -I/usr/local/include -DH5_USE_16_API -I./ 
INCL += -L/usr/local/lib -lhdf5 -lz #-ltbb -D_GLIBCXX_USE_TBB_PAR_BACKEND=0


all:
	@make falcon

falcon: main.cxx
	$(CXX) $? -o $@ $(OPT) $(INCL)
	#./Taichi 0 dat.10 500 1
	#./Taichi 1 snapshot_2.hdf5 500 1

clean:
	$(RM) ./*.o ./falcon ./fmm ./*~
