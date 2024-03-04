#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace std;
using namespace Rcpp;

using PermMat = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>;

// [[Rcpp::export]]
void test_permute_matrix()
{
    MatrixXd m(2, 5);
    m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
    cout << "m: before permutation" << endl;
    cout << m << endl;
    VectorXi indices(5);
    indices << 3, 0, 4, 2, 1;
    PermMat p(indices);
    cout << "p: permutation" << endl;
    cout << p.indices().transpose() << endl;
    cout << "m: after permutation by columns" << endl;
    m = m * p;
    cout << m << endl;
    cout << "m: recover the origianl by t(p*t(m))" << endl;
    m = (p * m.transpose()).transpose();
    cout << m << endl;
    cout << "o: the ordering for mapping the original to the permuted matrix" << endl;
    cout << "use p to get o" << endl;
    VectorXi ord(5);
    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++)
            if(indices[j] == i) {
                ord[i] = j;break;
            }
    }
    // ord << 1, 4, 3, 0, 2;
    cout << ord.transpose() << endl;
    cout << "the relationship between p and o" << endl;
    for(int i = 0; i < 5; i++) cout << "i:"<< i << ", p:" << indices[i] << ", o:" << ord[i] << "\n";
    cout << "p[o[i]] = i" << endl;
    for(int i = 0; i < 5; i++) cout << p.indices()[ord[i]] << "\n";
}

