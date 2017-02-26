#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <opencv2/core/eigen.hpp>

using namespace Eigen;
using namespace std;
using namespace cv;

//void fastLaplace(MatrixXd Dx, MatrixXd Dy, MatrixXd temp, MatrixXd& Lx, MatrixXd& Ly);



void fastLaplace(MatrixXd Dx, MatrixXd Dy, MatrixXd temp, MatrixXd& Lx, MatrixXd& Ly) {
	int a = Dx.cols();
	int b = Dx.rows();
	int n = a*b;
	int m = a*b;

	///temp == rotatedBoundary(:, : , 1)
	///Dx == rotatedBoundary(:, : , 2)
	///Dy == rotatedBoundary(:, : , 3)

	VectorXd AA(m);
	VectorXd BB(m);
	VectorXd i(m * 5);
	VectorXd j(m * 5);
	VectorXd s(m * 5);

	double t = 0;
	double k;

	for (int p = 0; p < m; p++) {
		double X = floor((p - 1) / a) + 1;
		double Y = p - (X - 1)*a;
		if (temp(Y, X) > 0) {
			i(t) = p;
			j(t) = p;
			s(t) = 1;
			t = t + 1;
			AA(p) = Dx(Y, X);
			BB(p) = Dy(Y, X);
		}
		else {
			i(t) = p;
			j(t) = p;
			s(t) = 1;
			t = t + 1;

			k = 0;
			if ((p - a) > 0)
				k = k + 2;
			else if ((p - 1) > 0)
				k = k + 1;

			if ((p + a) <= m)
				k = k + 2;
			else if ((p + 1) <= m)
				k = k + 1;

			k = double(k);


			if ((p - a) > 0) {
				i(t) = p;
				j(t) = p - a;
				s(t) = double(-1.0 / k);
				t = t + 1;
				i(t) = p;
				j(t) = p - 1;
				s(t) = double(-1.0 / k);
				t = t + 1;
			}
			else if ((p - 1) > 0) {
				i(t) = p;
				j(t) = p - 1;
				s(t) = double(-1.0 / k);
				t = t + 1;
			}
			if ((p + a) <= m) {
				i(t) = p;
				j(t) = p + 1;
				s(t) = double(-1.0 / k);
				t = t + 1;
				i(t) = p;
				j(t) = p + a;
				s(t) = double(-1.0 / k);
				t = t + 1;
			}
			else if ((p + 1) <= m) {
				i(t) = p;
				j(t) = p + 1;
				s(t) = double(-1.0 / k);
				t = t + 1;
			}
		}//end else
	}//end first for

	SparseMatrix<double> A(m, m);
	SparseMatrix<double> B(m, m);
	for (int ind1 = 0; ind1 < t-1; ind1++) {
		A.coeffRef(i(ind1), j(ind1)) = s(ind1);
		B.coeffRef(i(ind1), j(ind1)) = s(ind1);
	}
	//A = sparse(i(1:t - 1), j(1:t - 1), s(1:t - 1), m, m, m * 5);
	//B = sparse(i(1:t - 1), j(1:t - 1), s(1:t - 1), m, m, m * 5);

	// Solving:
	SimplicialCholesky<SparseMatrix<double>> cholA(A);  // performs a Cholesky factorization of A
	VectorXd LA = cholA.solve(AA);         // use the factorization to solve for the given right hand side

	SimplicialCholesky<SparseMatrix<double>> cholB(B);  // performs a Cholesky factorization of A
	VectorXd LB = cholB.solve(BB);         // use the factorization to solve for the given right hand side

	//LA = A\AA;
	//LB = B\BB;

	for (int e = 0; e < a; e++) {
		for (int w = 0; w < b; b++) {
			Ly(e, w) = LB(e + (w - 1)*a);
			Lx(e, w) = LA(e + (w - 1)*a);
		}
	}
}