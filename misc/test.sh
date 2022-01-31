#!/usr/bin/env bash
#

# download examples data

if [ `which curl` ];then
    curl -O http://popgen.dk/zilong/datahub/pca/examples.tar.gz
elif [ `which wget` ];then
    wget http://popgen.dk/zilong/datahub/pca/examples.tar.gz
else
    echo "please make sure curl or wget installed"
    exit 1
fi

tar -xzf examples.tar.gz

K=3

chmod +x PCAone

./PCAone --bfile examples/asia -v -k $K -o examples/test.arnoldi && echo "PCAone Arnoldi batch mode ok"

./PCAone --bfile examples/asia -v -k $K -o examples/test.arnoldi.m -m 1 && echo "PCAone Arnoldi block mode ok"

./PCAone --bfile examples/asia -v -k $K -o examples/test.halko -h && echo "PCAone Halko batch mode ok"

./PCAone --bfile examples/asia -v -k $K -o examples/test.halko.m -h -m 1 && echo "PCAone Halko block mode ok"

./PCAone --bfile examples/asia -v -k $K -o examples/test.fast.halko -f && echo "PCAone Fast Halko batch mode ok"

./PCAone --bfile examples/asia -v -k $K -o examples/test.fast.halko.m -f -m 1 && echo "PCAone Fast Halko block mode ok"

./PCAone --bgen examples/asia.bgen -v -k $K -o examples/test.bgen.arnoldi && echo "PCAone FileBgen batch mode ok"

./PCAone --beagle examples/sim.bgl.gz -p -v -k $K -o examples/test.bgl.arnoldi && echo "PCAone FileBeagle batch mode ok"

./PCAone --csv examples/data.csv.zst -v -k $K -o examples/test.csv.arnoldi && echo "PCAone FileCsv batch mode ok"
