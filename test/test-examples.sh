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

./PCAone --bfile examples/test -v -n 2 -k $K -o examples/test.arnoldi && echo "PCAone Arnoldi batch mode ok"

# ./PCAone --bfile examples/test -v -n 2 -k $K -o examples/test.arnoldi.m -m 1 && echo "PCAone Arnoldi block mode ok"

./PCAone --bfile examples/test -v -n 2 -k $K -o examples/test.halko -h --maxp 10 && echo "PCAone Halko batch mode ok"

# ./PCAone --bfile examples/test -v -n 2 -k $K -o examples/test.halko.m -h -m 1 && echo "PCAone Halko block mode ok"

./PCAone --bfile examples/test -v -n 2 -k $K -o examples/test.fast.halko -f --maxp 10 --shuffle && echo "PCAone Fast Halko batch mode ok"

# ./PCAone --bfile examples/test -v -n 2 -k $K -o examples/test.fast.halko.m -f -m 1 && echo "PCAone Fast Halko block mode ok"

# ./PCAone --bgen examples/test.bgen -v -n 2 -k $K -o examples/test.bgen.arnoldi && echo "PCAone FileBgen batch mode ok"

./PCAone --beagle examples/test.bgl.gz -p -n 2 -v -k $K -o examples/test.bgl.arnoldi && echo "PCAone FileBeagle batch mode ok"

# ./PCAone --csv examples/data.csv.zst -v -n 2 -k $K -o examples/test.csv.arnoldi && echo "PCAone FileCsv batch mode ok"

# ./PCAone --csv examples/data.csv.zst -v -n 2 -k $K -o examples/test.csv.arnoldi -m 0.01 && echo "PCAone FileCsv out-of-core mode ok"
