#!/usr/bin/env bash
#

# download examples data

if [ `which curl` ];then
    curl -O http://popgen.dk/zilong/datahub/pca/examples.tar.gz
elif [ `which wget` ];then
    wget -O examples.tar.gz http://popgen.dk/zilong/datahub/pca/examples.tar.gz
else
    echo "please make sure curl or wget installed"
    exit 1
fi

tar -xzf examples.tar.gz

K=3

if [ ! -x PCAone ]; then
    curl -LO https://github.com/Zilong-Li/PCAone/releases/latest/download/PCAone-avx2-Linux.zip && unzip PCAone-avx2-Linux.zip
fi

chmod +x PCAone

./PCAone --bfile examples/asia -f -v -k 3 -o out && cat out.eigvals && echo "PCAone fancy batch mode ok"

./PCAone --bfile examples/asia -f -v -k 3 -o out --printv && cat out.eigvals && echo "PCAone fancy batch mode --printv ok"

./PCAone --bfile examples/asia -f -v -k 3 -o out -m 2 && cat out.eigvals && echo "PCAone fancy block mode ok"

./PCAone --bfile examples/asia -f -v -k 3 -o out -m 2 --printv && cat out.eigvals && echo "PCAone fancy block mode -- printv ok"

./PCAone --bfile examples/asia -a -v -k 3 -o out && cat out.eigvals && echo "PCAone Arnoldi batch mode ok"

./PCAone --bfile examples/asia -a -v -k 3 -o out -m 2 && cat out.eigvals && echo "PCAone Arnoldi block mode ok"

./PCAone --bfile examples/asia -h -v -k 3 -o out && cat out.eigvals && echo "PCAone RSVD batch mode ok"

./PCAone --bfile examples/asia -h -v -k 3 -o out -m 2 && cat out.eigvals && echo "PCAone RSVD block mode ok"

./PCAone --csv examples/test.csv.zst -k 20 -n 10 -o out -a --cpmed && cat out.eigvals && echo "PCAone CSV batch mode ok"

./PCAone --bfile examples/test.emu -k 3 -n 10 -o out -m 2 --emu && cat out.eigvals && echo "PCAone EMU mode ok"

./PCAone --beagle examples/test.bgl.gz -k 3 -n 10 -o out -a --pcangsd && cat out.eigvals && echo "PCAone PCAngsd batch mode ok"
