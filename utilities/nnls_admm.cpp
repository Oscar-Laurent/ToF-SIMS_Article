#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Eigen/Dense>

namespace py = pybind11;
using namespace Eigen;

std::pair<py::array_t<double>, double> nnls_admm(py::array_t<double> A_in, py::array_t<double> b_in, py::array_t<double> Y_in, 
                                                 double rho = 1.0, double alpha = 1.0, int max_iter = 1000, double tol = 1e-8) {
    // Convert input arrays to Eigen matrices
    Eigen::Map<MatrixXd> A(A_in.mutable_data(), A_in.shape(0), A_in.shape(1));
    Eigen::Map<VectorXd> b(b_in.mutable_data(), b_in.size());
    Eigen::Map<VectorXd> Y(Y_in.mutable_data(), Y_in.size());

    int m = A.rows();
    int n = A.cols();

    VectorXd x = VectorXd::Zero(n);
    VectorXd z = VectorXd::Zero(m);
    VectorXd u = VectorXd::Zero(m);

    MatrixXd AtA = A.transpose() * A;
    MatrixXd I = MatrixXd::Identity(n, n);

    Eigen::LLT<MatrixXd> llt(AtA + rho * I);

    for (int k = 0; k < max_iter; ++k) {
        // x-update
        VectorXd q = A.transpose() * (b + z - u) + rho * x;
        x = llt.solve(q);
        x = x.cwiseMax(0.0);

        // z-update with over-relaxation
        VectorXd Ax = A * x;
        VectorXd z_old = z;
        z = (Ax + u).cwiseMin(Y);
        z = alpha * z + (1 - alpha) * z_old;

        // u-update
        u += Ax - z;

        // Check convergence
        if ((z - z_old).norm() < tol) {
            break;
        }
    }

    double rnorm = (A * x - b).norm();
    
    return {py::array_t<double>(x.size(), x.data()), rnorm};
}

PYBIND11_MODULE(nnls_admm, m) {
    m.def("nnls_admm", &nnls_admm, "Solve NNLS with constraints using ADMM",
          py::arg("A"), py::arg("b"), py::arg("Y"),
          py::arg("rho") = 1.0, py::arg("alpha") = 1.0, py::arg("max_iter") = 1000, py::arg("tol") = 1e-8);
}