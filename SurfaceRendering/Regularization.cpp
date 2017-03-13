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
#include <Eigen/Sparse>
#include <stdio.h>      /* printf */
#include <math.h>  
#include <unsupported/Eigen/NonLinearOptimization>
//#include <unsupported/Eigen/LevenbergMarquardt>


// LM minimize for the model y = a x + b
typedef std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > Point2DVector;

using namespace Eigen;
using namespace std;
using namespace cv;

int overlap = 1;
int counter = 0;

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
	typedef _Scalar Scalar;
	enum {
		InputsAtCompileTime = NX,
		ValuesAtCompileTime = NY
	};
	typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
	typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
	typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

	int m_inputs, m_values;

	Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
	Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

	int inputs() const { return m_inputs; }
	int values() const { return m_values; }

};


struct MyFunctor : Functor<double>
{
	int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
	{
		// "a" in the model is x(0), and "b" is x(1)
		for (unsigned int i = 0; i < this->Points.size(); ++i)
		{
			//fvec(i) = this->Points[i](1) - (x(0) * sin(x(1)*this->Points[i](0) + x(2)) + x(3)); for sine-fitting
			fvec(i) = this->Points[i](1) - (x(0)); //for constant-fitting
		}

		return 0;
	}

	Point2DVector Points;

	int inputs() const { return 1; } // There are two parameters of the model
	int values() const { return this->Points.size(); } // The number of observations
};

struct MyFunctorNumericalDiff : Eigen::NumericalDiff<MyFunctor> {};

Point2DVector GeneratePoints(const unsigned int numberOfPoints)
{
	Point2DVector points;
	// Model y = 2*x + 5 with some noise (meaning that the resulting minimization should be about (2,5)
	for (unsigned int i = 0; i < numberOfPoints; ++i)
	{
		double x = static_cast<double>(i);
		Eigen::Vector2d point;
		point(0) = x;
		point(1) = 50 * sin(2 * x + 1) + 50;
		points.push_back(point);
	}

	return points;
}


Mat warpRegression(Mat im, int startCol, int endCol, int startRow, int endRow) {

	MatrixXd im_mat;
	cv2eigen(im, im_mat);

	Point2DVector points;
	int indx, col =0;
	//TO DO: make sure you changed the size of image back
	for (int i = startCol; i <= endCol; i++) {
		for (int j = 0; j < im_mat.rows(); j++) {
			indx = j + col*im_mat.rows();
			Eigen::Vector2d point;
			point(0) = indx;
			point(1) = im_mat(j,i);
			points.push_back(point);

		}
		col++;
	}

	//unsigned int numberOfPoints = 200;
	//Point2DVector points = GeneratePoints(numberOfPoints);

	//initialize the theta vector
	Eigen::VectorXd theta(1);
	theta << 233 / 2;
	//x.fill(4.0f);

	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
	//std::cout << "status: " << status << std::endl;

	MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

	int min_col, max_col, min_row, max_row;
	if (startCol == 0 && startRow == 0) {
		//cout << "case1" << endl;
		min_col = startCol;
		min_row = startRow;
		max_col = endCol + overlap;
		max_row = endRow - overlap;
	}
	else if (endCol == 456 && startRow == 0) {
		//cout << "case2" << endl;
		max_col = endCol;
		min_row = startRow;
		min_col = startCol - overlap;
		max_row = endRow + overlap;
	}
	else if (startCol == 0 && endRow == 670) {
		//cout << "case3" << endl;
		min_col = startCol;
		max_row = endRow;
		max_col = endCol + overlap;
		min_row = startRow - overlap;
	}
	else if (endCol == 456 && endRow == 670) {
		//cout << "case 4" << endl;
		max_col = endCol;
		max_row = endRow;
		min_col = startCol - overlap;
		min_row = startRow - overlap;
	}
	else {
		min_col = startCol - overlap;
		min_row = startRow - overlap;
		max_col = endCol + overlap;
		max_row = endRow + overlap;
	}

	for (int i = min_col  ; i < max_col ; i++) {
		for (int j = min_row ; j < max_row ; j++) {
			double x = j;
			double y = theta(0);
			im_reg_mat(j, i) = y;
		}
	}


	Mat regIm, temp;
	cv::eigen2cv(im_reg_mat, temp);
	temp.convertTo(regIm, CV_8UC1);
	//imshow("regularized", regIm);
	cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
	std::cout << "theta for patch_" + std::to_string(counter) << ":  " << theta << std::endl;
	counter++;

	return regIm;
}

Mat weftRegression(Mat im, int startCol, int endCol, int startRow, int endRow) {

	MatrixXd im_mat;
	cv2eigen(im, im_mat);

	Point2DVector points;
	int indx, row = 0;
	//TO DO: make sure you changed the size of image back

	for (int i = startRow; i <= endRow; i++) {
		for (int j = startCol; j < endCol; j++) {
			indx = j + row*im_mat.cols();
			Eigen::Vector2d point;
			point(0) = indx;
			point(1) = im_mat(i,j);
			points.push_back(point);

		}
		row++;
	}


	//initialize the theta vector
	Eigen::VectorXd theta(1);
	theta << 233 / 3;
	//x.fill(4.0f);

	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
	//std::cout << "status: " << status << std::endl;

	MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

	int min_col, max_col, min_row, max_row;
	if (startCol == 1 && startRow == 1) {
		//cout << "case1" << endl;
		min_col = startCol;
		min_row = startRow;
		max_col = endCol + overlap;
		max_row = endRow - overlap;
	}
	else if (endCol == 456 && startRow == 1) {
		//cout << "case2" << endl;
		max_col = endCol;
		min_row = startRow;
		min_col = startCol - overlap;
		max_row = endRow + overlap;
	}
	else if (startCol == 1 && endRow == 670) {
		//cout << "case3" << endl;
		min_col = startCol;
		max_row = endRow;
		max_col = endCol + overlap;
		min_row = startRow - overlap;
	}
	else if (endCol == 456 && endRow == 670) {
		//cout << "case 4" << endl;
		max_col = endCol;
		max_row = endRow;
		min_col = startCol - overlap;
		min_row = startRow - overlap;
	}
	else {
		//cout << "case 5" << endl;
		min_col = startCol - overlap;
		min_row = startRow - overlap;
		max_col = endCol + overlap;
		max_row = endRow + overlap;
	}


	for (int i = min_row; i < max_row ; i++) {
		for (int j = min_col ; j < max_col ; j++) {
			double x = j;
			double y = theta(0);
			im_reg_mat(i, j) = y;
		}
	}

	Mat regIm, temp;
	cv::eigen2cv(im_reg_mat, temp);
	temp.convertTo(regIm, CV_8UC1);
	//imshow("regularized", regIm);
	cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
	std::cout << "theta for patch_" + std::to_string(counter) << ":  " << theta << std::endl;
	counter++;

	return regIm;
}

vector<Mat> regularization(vector<Mat> morphed_patches, vector < vector<Point> > fixedPoints, int padding);
vector<Mat> regularization(vector < vector<Point> > fixedPoints, int padding ) {
	vector<Mat> morphed_patches;
	for (int i = 0; i < 30; i++) {
		morphed_patches.push_back(imread("Morphed Patches/morphed_patch_" + std::to_string(i) + ".png", CV_LOAD_IMAGE_GRAYSCALE));
	}
	vector<Mat> reg_morphed_patches = regularization(morphed_patches, fixedPoints, padding);
	return reg_morphed_patches;
}
vector<Mat> regularization(vector<Mat> morphed_patches, vector < vector<Point> > fixedPoints, int padding) {

	vector<Mat> reg_morphed_patches;

	for (int i = 0; i < 30; i++) {

		reg_morphed_patches.push_back(weftRegression(morphed_patches[i], fixedPoints[i][0].x - padding, fixedPoints[i][3].x - padding, fixedPoints[i][0].y - padding, fixedPoints[i][3].y - padding));
		i++;
		reg_morphed_patches.push_back(warpRegression(morphed_patches[i], fixedPoints[i][0].x - padding, fixedPoints[i][3].x - padding, fixedPoints[i][0].y - padding, fixedPoints[i][3].y - padding) );	
		
	}

	return reg_morphed_patches;
}