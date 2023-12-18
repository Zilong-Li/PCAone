#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
VectorXd rcppeigen_cor(MatrixXd  X) {
    X.rowwise() -= X.colwise().mean(); // Centering
    const int df = X.rows() - 1; // N-1
    VectorXd sds = (X.array().square().colwise().sum() / df).sqrt();
    VectorXd res(X.cols());
    for(Eigen::Index i = 0; i < X.cols(); i++) {
        res(i)  = (X.col(0).array() * X.col(i).array() / (sds(0) * sds(i))).sum() / df;
    }
    return res;
}

// [[Rcpp::export]]
vector<int> prune(MatrixXd X, VectorXi pos, int window, double cutoff) {
    X.rowwise() -= X.colwise().mean(); // Centering
    const int df = X.rows() - 1; // N-1
    VectorXd sds = (X.array().square().colwise().sum() / df).sqrt();
    int m = X.cols();
    vector<int> keep(m, 1);
    Rcout << "test prune\n";
    for(int i = 0; i < m; i++) {
        if(!keep[i]) continue;
        for(int j = i + 1; j < m; j++) {
            if(!keep[j]) continue;
            if(pos[j] - pos[i] > window) break;
            double r  = (X.col(j).array() * X.col(i).array() / (sds(j) * sds(i))).sum() / df;
            if( r * r > cutoff) keep[j] = 0;
        }
    }
    
    return keep;
}
