#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <opencv2/core/eigen.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include "ImageProcessing.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

using namespace Eigen;
using namespace std;
using namespace cv;
typedef Eigen::SparseMatrix <double> SpMat;


MatrixXd fastLaplace(MatrixXd Dx, MatrixXd cp) {
	int height = Dx.rows();
	int width = Dx.cols();
	MatrixXd Dx_pad = MatrixXd::Zero(height + 2, width + 2);
	int height_pad = Dx_pad.rows();
	int width_pad = Dx_pad.cols();
	//MatrixXd A = MatrixXd::Identity(height_pad*width_pad, height_pad*width_pad);
	SpMat A(height_pad*width_pad, height_pad*width_pad);
	A.reserve(5 * height_pad*width_pad);


	MatrixXd id(height_pad, width_pad);
	VectorXd b = VectorXd::Zero(height_pad*width_pad);
	MatrixXd Lx_pad(height_pad, width_pad);
	MatrixXd Lx(height, width);

	///Zero padding the input
	for (int i = 1; i <= height; i++) {
		for (int j = 1; j <= width; j++) {
			Dx_pad(i, j) = Dx(i - 1, j - 1);
		}
	}

	///define id matrix to convert Dx_pad and Dy_pad to a vector
	int n = 0;
	for (int i = 0; i < height_pad; i++) {
		for (int j = 0; j < width_pad; j++) {
			id(i, j) = n;
			n++;
		}
	}

	/// Create matrix A and Vector b
	for (int i = 1; i <= height; i++) {
		for (int j = 1; j <= width; j++) {
			if (cp(i - 1, j - 1) < 1) {
				int top_id = id(i - 1, j);
				int bottom_id = id(i + 1, j);
				int right_id = id(i, j + 1);
				int left_id = id(i, j - 1);
				cout << top_id << "  " << bottom_id << "  " << right_id << "  " << left_id << endl;
				A.insert(id(i, j), top_id) = -0.25;
				A.insert(id(i, j), bottom_id) = -0.25;
				A.insert(id(i, j), right_id) = -0.25;
				A.insert(id(i, j), left_id) = -0.25;
			}
			b(id(i, j)) = Dx_pad(i, j);
		}
	}

	///solve the linear system Ax = b
	A.makeCompressed();
	SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solverA;
	solverA.compute(A);
	VectorXd x = solverA.solve(b);
	///alternative ways:
	//VectorXd x = A.colPivHouseholderQr().solve(b);
	// or
	//Eigen::SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
	//Eigen::VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side

	///convert vector x to matrix
	n = 0;
	for (int i = 0; i < height_pad; i++) {
		for (int j = 0; j < width_pad; j++) {
			Lx_pad(i, j) = x(n);
			n++;
		}
	}

	///crop the margin
	for (int i = 1; i <= height; i++) {
		for (int j = 1; j <= width; j++) {
			Lx(i - 1, j - 1) = Lx_pad(i, j);
		}
	}

	return Lx;

#if 0
	///check how to solve a sparce matrix linear system:
	SparseMatrix<double> A(3, 3);
	A.insert(0, 0) = 1;
	A.insert(0, 1) = 1;
	A.insert(0, 2) = 1;
	A.insert(1, 0) = 0;
	A.insert(1, 1) = 2;
	A.insert(1, 2) = 5;
	A.insert(2, 0) = 2;
	A.insert(2, 1) = 5;
	A.insert(2, 2) = -1;

	A.makeCompressed();
	SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solverA;
	solverA.compute(A);

	VectorXd B;
	B.resize(3);
	B << 6, -4, 27;

	VectorXd X = solverA.solve(B);

	cout << A << endl << B << endl << X << endl;
#endif
}

/// solving Laplace's equation using the Jacobi method.
void slowLaplace(MatrixXd Dx, MatrixXd Dy, MatrixXd temp, MatrixXd& Lx, MatrixXd& Ly) {
	int n = Dx.rows();
	int m = Dx.cols();
	float tol = 1000;
	float err = 1000;
	int k = 0;

	Lx = Dx;
	///iterate Jacobi until convergence
	while (err > tol) {
		k++;
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < m - 1; j++) {
				if (temp(i, j)) {
					Lx(i, j) = Dx(i, j);
				}
				else {
					Lx(i, j) = (Dx(i - 1, j) + Dx(i + 1, j) + Dx(i, j - 1) + Dx(i, j + 1))*0.25;
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
	
	
	slowLaplace(dx, dy, temp, Lx, Ly);
	//Lx = fastLaplace(dx, temp);
	//Ly = fastLaplace(dy, temp);

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
	cout << "morphing() is returned ... " << endl;
	eigen2cv(img_new_mat, img_new);

	return;
}
