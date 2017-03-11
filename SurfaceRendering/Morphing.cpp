#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <opencv2/core/eigen.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include "ImageProcessing.h"

using namespace Eigen;
using namespace std;
using namespace cv;


//fastLaplace(Dx, Dy, temp, Lx, Ly);

/// solving Laplace's equation using the Jacobi method.
void fastLaplace(MatrixXd Dx, MatrixXd Dy, MatrixXd temp, MatrixXd& Lx, MatrixXd& Ly) {
	int n = Dx.rows();
	int m = Dx.cols();
	float tol = 1000;
	float err = 1000;
	int k = 0;

	Lx = Dx;
	///iterate Jacobi until convergence
	while (err > tol ) {
		k++;
		for (int i = 1; i < n - 1; i++){
			for (int j = 1; j < m - 1; j++) {
				if (temp(i, j)) {
					Lx(i, j) = Dx(i, j);
				}
				else {
					Lx(i, j) = (Dx(i - 1, j) + Dx(i + 1, j) + Dx(i, j - 1) + Dx(i, j + 1) )*0.25;
				}
			}
		}
		err = ((Lx - Dx).cwiseAbs()).sum(); //L1 norm
		Dx = Lx;
	}

	k = 0;
	err = 1000;
	Ly = Dy;
	///iterate Jacobi until convergence
	while (err > tol) {
		k++;
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < m - 1; j++) {
				if (temp(i, j)) {
					Ly(i, j) = Dy(i, j);
				}
				else {
					Ly(i, j) = (Dy(i - 1, j) + Dy(i + 1, j) + Dy(i, j - 1) + Dy(i, j + 1))*0.25;
				}
			}
		}
		err = ((Ly - Dy).cwiseAbs()).sum(); //L1 norm
		Dy = Ly;
	}

	return;
}

void morphing(Mat& img_new, Mat img, vector<Point> movingPoints, vector<Point> fixedPoints) {

	MatrixXd img_mat;
	cv2eigen(img, img_mat);
	int m = img_mat.rows();
	int n = img_mat.cols();
	MatrixXd img_new_mat = MatrixXd::Zero(m, n);
	MatrixXd dx = MatrixXd::Zero(m, n);
	MatrixXd dy = MatrixXd::Zero(m, n);
	MatrixXd temp = MatrixXd::Zero(m, n);
	int s = movingPoints.size();

	for (int i = 0; i < s; i++) {
		dx(fixedPoints[i].y, fixedPoints[i].x) = movingPoints[i].x - fixedPoints[i].x;
		dy(fixedPoints[i].y, fixedPoints[i].x) = movingPoints[i].y - fixedPoints[i].y;
		temp(fixedPoints[i].y, fixedPoints[i].x) = 1;
	}

	MatrixXd Lx, Ly;
	fastLaplace(dx, dy, temp, Lx, Ly);

	///display the Lx and Ly
	//Mat t;
	//eigen2cv(Lx, t);
	//imshow("Lx", t);
	//eigen2cv(Ly, t);
	//imshow("Ly", t);

	int padding = 100;
	Mat img_pad = ZeroPadding(img, padding);
	MatrixXd img_pad_mat;
	cv2eigen(img_pad, img_pad_mat);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			img_new_mat(i, j) = img_pad_mat(i + Lx(i, j) + padding, j + Ly(i, j) + padding);
		}
	}
	//TO DO change the size back
	cout <<  "morphing() is returned ... " << endl;
	eigen2cv(img_new_mat, img_new);

	return;
}

//void fastLaplace(MatrixXd Dx, MatrixXd Dy, MatrixXd temp, MatrixXd& Lx, MatrixXd& Ly) {
//	int a = Dx.cols();
//	int b = Dx.rows();
//	int n = a*b;
//	int m = a*b;
//
//	///temp == rotatedBoundary(:, : , 1)
//	///Dx == rotatedBoundary(:, : , 2)
//	///Dy == rotatedBoundary(:, : , 3)
//
//	VectorXd AA(m);
//	VectorXd BB(m);
//	VectorXd i(m * 5);
//	VectorXd j(m * 5);
//	VectorXd s(m * 5);
//
//	double t = 0;
//	double k;
//
//	for (int p = 0; p < m; p++) {
//		double X = floor((p - 1) / a) + 1;
//		double Y = p - (X - 1)*a;
//		if (temp(Y, X) > 0) {
//			i(t) = p;
//			j(t) = p;
//			s(t) = 1;
//			t = t + 1;
//			AA(p) = Dx(Y, X);
//			BB(p) = Dy(Y, X);
//		}
//		else {
//			i(t) = p;
//			j(t) = p;
//			s(t) = 1;
//			t = t + 1;
//
//			k = 0;
//			if ((p - a) > 0)
//				k = k + 2;
//			else if ((p - 1) > 0)
//				k = k + 1;
//
//			if ((p + a) <= m)
//				k = k + 2;
//			else if ((p + 1) <= m)
//				k = k + 1;
//
//			k = double(k);
//
//
//			if ((p - a) > 0) {
//				i(t) = p;
//				j(t) = p - a;
//				s(t) = double(-1.0 / k);
//				t = t + 1;
//				i(t) = p;
//				j(t) = p - 1;
//				s(t) = double(-1.0 / k);
//				t = t + 1;
//			}
//			else if ((p - 1) > 0) {
//				i(t) = p;
//				j(t) = p - 1;
//				s(t) = double(-1.0 / k);
//				t = t + 1;
//			}
//			if ((p + a) <= m) {
//				i(t) = p;
//				j(t) = p + 1;
//				s(t) = double(-1.0 / k);
//				t = t + 1;
//				i(t) = p;
//				j(t) = p + a;
//				s(t) = double(-1.0 / k);
//				t = t + 1;
//			}
//			else if ((p + 1) <= m) {
//				i(t) = p;
//				j(t) = p + 1;
//				s(t) = double(-1.0 / k);
//				t = t + 1;
//			}
//		}//end else
//	}//end first for
//
//	SparseMatrix<double> A(m, m);
//	SparseMatrix<double> B(m, m);
//	for (int ind1 = 0; ind1 < t-1; ind1++) {
//		A.coeffRef(i(ind1), j(ind1)) = s(ind1);
//		B.coeffRef(i(ind1), j(ind1)) = s(ind1);
//	}
//	//A = sparse(i(1:t - 1), j(1:t - 1), s(1:t - 1), m, m, m * 5);
//	//B = sparse(i(1:t - 1), j(1:t - 1), s(1:t - 1), m, m, m * 5);
//
//	// Solving:
//	SimplicialCholesky<SparseMatrix<double>> cholA(A);  // performs a Cholesky factorization of A
//	VectorXd LA = cholA.solve(AA);         // use the factorization to solve for the given right hand side
//
//	SimplicialCholesky<SparseMatrix<double>> cholB(B);  // performs a Cholesky factorization of A
//	VectorXd LB = cholB.solve(BB);         // use the factorization to solve for the given right hand side
//
//	//LA = A\AA;
//	//LB = B\BB;
//
//	for (int e = 0; e < a; e++) {
//		for (int w = 0; w < b; b++) {
//			Ly(e, w) = LB(e + (w - 1)*a);
//			Lx(e, w) = LA(e + (w - 1)*a);
//		}
//	}
//}