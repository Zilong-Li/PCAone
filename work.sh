g++ -O3 -Wall -march=native -std=c++11 -DEIGEN_USE_BLAS -o emu EMU.cpp Data.cpp SVD.cpp Utils.cpp -I/home/zilong/local/include -I/usr/include/x86_64-linux-gnu -I./ -lopenblas
# g++ -O3 -Wall -march=native -std=c++11 -funroll-loops -ftree-vectorize -ffast-math -o emu EMU.cpp Data.cpp SVD.cpp -I/home/zilong/local/include -I/usr/include/x86_64-linux-gnu -I./ -lopenblas -lpthread -lgfortran
./emu /home/jonas/bgi/siyang/newSimulation/sim.nLow.mMed.geno 10
