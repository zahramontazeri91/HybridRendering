#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/MatrixFunctions>
#include "ImageProcessing.h"
#include "input.h"
#include <Eigen/Sparse>
#include <stdio.h>      /* printf */
#include <math.h>  
#include <unsupported/Eigen/NonLinearOptimization>
//#include <unsupported/Eigen/LevenbergMarquardt>
#include <opencv2/core/eigen.hpp>

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
			//fvec(i) = this->Points[i](1) - (x(0)); //for constant-fitting 
			fvec(i) = this->Points[i](1) - (x(0) * tanh(x(1)*(this->Points[i](0) + x(2))) + x(3)); //two back-to-back tanh fitting
		}

		return 0;
	}

	Point2DVector Points;

	int inputs() const { return 4; } // There are two parameters of the model
	int values() const { return this->Points.size(); } // The number of observations
};

struct MyFunctorNumericalDiff : Eigen::NumericalDiff<MyFunctor> {};


Mat warpRegression(Mat im, int startCol, int endCol, int startRow, int endRow) {


	int min_col, max_col, min_row, max_row;
	if (startCol == 1) {
		min_col = startCol;
		min_row = startRow - overlap;
		max_col = endCol + overlap;
		max_row = endRow + overlap;

		if (startRow == 1) {
			min_row = startRow;
		}
		if (endRow == 670) {
			max_row = endRow;
		}
	}
	else if (endCol == 456) {
		max_col = endCol;
		min_row = startRow - overlap;
		min_col = startCol - overlap;
		max_row = endRow + overlap;

		if (startRow == 1) {
			min_row = startRow;
		}
		if (endRow == 670) {
			max_row = endRow;
		}
	}
	else if (startRow == 1) {
		min_row = startRow;
		min_col = startCol - overlap;
		max_col = endCol + overlap;
		max_row = endRow + overlap;

		if (startCol == 1) {
			min_col = startCol;
		}
		if (endCol == 456) {
			max_col = endCol;
		}
	}
	else if (endRow == 670) {
		max_row = endRow;
		min_row = startRow - overlap;
		min_col = startCol - overlap;
		max_col = endCol + overlap;

		if (startCol == 1) {
			min_col = startCol;
		}
		if (endCol == 456) {
			max_col = endCol;
		}
	}
	else {
		min_col = startCol - overlap;
		min_row = startRow - overlap;
		max_col = endCol + overlap;
		max_row = endRow + overlap;
	}


	MatrixXd im_mat;
	cv2eigen(im, im_mat);

	Point2DVector points;
	int indx, col =0;

	int transition=  int( (30.0 / 100.0) * (endRow - startRow) );
 
	for (int i = startCol; i <= endCol; i++) {
		indx = 0;
		for (int j = startRow; j <= startRow + transition; j++) {
			//indx = j + col*(endRow - startRow) + 1;
			Eigen::Vector2d point;
			point(0) = indx;
			point(1) = im_mat(j,i);
			points.push_back(point);
			indx++;
		}
		col++;
	}

	//unsigned int numberOfPoints = 200;
	//Point2DVector points = GeneratePoints(numberOfPoints);

	//initialize the theta vector
	Eigen::VectorXd theta(4);
	double d = startRow + transition - min_row;
	//(theta(0)/2.0) * tanh(4.0/d* (x-d/2) ) + (theta(0)/2.0)
	theta << 170.0 / 2.0, 4.0/d, -1.0 * (d/2.0), 170.0 / 2.0;
	//x.fill(4.0f);

	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
	//std::cout << "status: " << status << std::endl;

	MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

	//cout << "min_row " << min_row << "startRow " << startRow << "startRow + transition " << startRow + transition << "endRow - transition " << endRow - transition << "endRow " << endRow << "max_row " << max_row << endl;

	double x;
	for (int i = min_col; i < max_col; i++) {
		x = 0;
		for (int j = min_row; j < startRow + transition; j++) {
			im_reg_mat(j, i) = (theta(0)) * tanh(theta(1)* (x + theta(2)) ) + (theta(3));
			if (j >= startRow) x++;   // in order to copy the first value for the overlapping region 
		}
	}

	// in order to copy the last value for the constant region in the middle 
	for (int i = min_col; i < max_col; i++) {
		for (int j = startRow + transition; j < endRow - transition; j++) {
			im_reg_mat(j, i) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
		}
	}

	for (int i = min_col; i < max_col; i++) {
		x = 0;
		for (int j = endRow - transition; j < max_row ; j++) {
			im_reg_mat(j, i) = (theta(0)) * tanh(-1 * theta(1)* (x + theta(2))) + (theta(3));
			if (j <= endRow) x++;
		}
	}

	Mat regIm, temp;
	cv::eigen2cv(im_reg_mat, temp);
	temp.convertTo(regIm, CV_8UC1);
	//imshow("regularized", regIm);
	cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
	std::cout << "theta for warp patch_" + std::to_string(counter) << ":  " << endl<< theta << std::endl;
	counter++;

	return regIm;
}

Mat weftRegression(Mat im, int startCol, int endCol, int startRow, int endRow) {

	int min_col, max_col, min_row, max_row;
	if (startCol == 1) {
		min_col = startCol;
		min_row = startRow - overlap;
		max_col = endCol + overlap;
		max_row = endRow + overlap;

		if (startRow == 1) {
			min_row = startRow;
		}
		if (endRow == 670) {
			max_row = endRow;
		}
	}
	else if (endCol == 456) {
		max_col = endCol;
		min_row = startRow - overlap;
		min_col = startCol - overlap;
		max_row = endRow + overlap;

		if (startRow == 1) {
			min_row = startRow;
		}
		if (endRow == 670) {
			max_row = endRow;
		}
	}
	else if (startRow == 1) {
		min_row = startRow;
		min_col = startCol - overlap;
		max_col = endCol + overlap;
		max_row = endRow + overlap;

		if (startCol == 1) {
			min_col = startCol;
		}
		if (endCol == 456) {
			max_col = endCol;
		}
	}
	else if (endRow == 670) {
		max_row = endRow;
		min_row = startRow - overlap;
		min_col = startCol - overlap;
		max_col = endCol + overlap;

		if (startCol == 1) {
			min_col = startCol;
		}
		if (endCol == 456) {
			max_col = endCol;
		}
	}
	else {
		min_col = startCol - overlap;
		min_row = startRow - overlap;
		max_col = endCol + overlap;
		max_row = endRow + overlap;
	}

	MatrixXd im_mat;
	cv2eigen(im, im_mat);

	Point2DVector points;
	int indx, row = 0;

	int transition = (20.0 / 100.0) * (endCol - startCol);
	for (int i = startRow; i <= endRow; i++) {
		indx = 0;
		for (int j = startCol; j <= startCol + transition; j++ ) {
			//indx = j + row* (endCol - startCol) + 1;
			Eigen::Vector2d point;
			point(0) = indx;
			point(1) = im_mat(i,j);
			points.push_back(point);
			indx++;
		}
		row++;
	}


	//initialize the theta vector
	Eigen::VectorXd theta(4);
	theta << 170 / 2, 4.0 / 44, -1.0* (44 / 2), 170 / 2;
	//x.fill(4.0f);

	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
	//std::cout << "status: " << status << std::endl;

	MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

	double d = startCol + transition - min_col;
	double x;
	for (int i = min_row; i < max_row; i++) {
		x = 0;
		for (int j = min_col; j < startCol + transition; j++) {
			//im_reg_mat(i, j) = (theta(0) / 2.0) * tanh( 4.0 / d* (x - d / 2)) + (theta(0) / 2.0);
			im_reg_mat(i, j) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
			if (j >= startCol) x++;	// in order to copy the first value for the overlapping region 
		}
	}

	// in order to copy the last value for the constant region in the middle 
	for (int i = min_row; i < max_row; i++) {
		for (int j = startCol + transition; j < endCol - transition; j++) {
			im_reg_mat(i, j) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
		}
	}

	for (int i = min_row; i < max_row ; i++) {
		x = 0;
		for (int j = endCol - transition ; j < max_col ; j++) {
			im_reg_mat(i, j) = (theta(0)) * tanh(-1 * theta(1)* (x + theta(2))) + (theta(3));
			if (j <=endCol) x++;	// in order to copy the first value for the overlapping region 
		}
	}

	Mat regIm, temp;
	cv::eigen2cv(im_reg_mat, temp);
	temp.convertTo(regIm, CV_8UC1);
	//imshow("regularized", regIm);
	cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
	std::cout << "theta for weft patch_" + std::to_string(counter) << ":  "<< endl << theta << std::endl;
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

	MatrixXd pattern = input_pattern();
	VectorXd columns = input_columns();
	VectorXd rows = input_rows();


	vector<Mat> reg_morphed_patches;

	for (int i = 0; i <= 3; i++) {

		reg_morphed_patches.push_back(weftRegression(morphed_patches[i], fixedPoints[i][0].x - padding, fixedPoints[i][3].x - padding, fixedPoints[i][0].y - padding, fixedPoints[i][3].y - padding));

		i++;

		if (i > 3) break;
		reg_morphed_patches.push_back(warpRegression(morphed_patches[i], fixedPoints[i][0].x - padding, fixedPoints[i][3].x - padding, fixedPoints[i][0].y - padding, fixedPoints[i][3].y - padding));
		

	}
	//cout << "first loop" << endl;
	for (int i = 4; i <= 16; i++) {

		reg_morphed_patches.push_back(warpRegression(morphed_patches[i], fixedPoints[i][0].x - padding, fixedPoints[i][3].x - padding, fixedPoints[i][0].y - padding, fixedPoints[i][3].y - padding));
		//cout << " p0 " << fixedPoints[i][0].x - padding << " " << fixedPoints[i][0].y - padding << endl <<
		//	" p1 " << fixedPoints[i][1].x - padding << "  " << fixedPoints[i][1].y - padding << endl <<
		//	" p2 " << fixedPoints[i][2].x - padding << "  " << fixedPoints[i][2].y - padding << endl <<
		//	" p3 " << fixedPoints[i][3].x - padding << "  " << fixedPoints[i][3].y - padding << endl;

		i++;

		if (i > 16) break;
		reg_morphed_patches.push_back(weftRegression(morphed_patches[i], fixedPoints[i][0].x - padding, fixedPoints[i][3].x - padding, fixedPoints[i][0].y - padding, fixedPoints[i][3].y - padding) );
		
	}
	//cout << "second loop" << endl;
	for (int i = 17 ; i <= 29; i++) {

		reg_morphed_patches.push_back(warpRegression(morphed_patches[i], fixedPoints[i][0].x - padding, fixedPoints[i][3].x - padding, fixedPoints[i][0].y - padding, fixedPoints[i][3].y - padding));

		i++;

		if (i > 29) break;
		reg_morphed_patches.push_back(weftRegression(morphed_patches[i], fixedPoints[i][0].x - padding, fixedPoints[i][3].x - padding, fixedPoints[i][0].y - padding, fixedPoints[i][3].y - padding));

	}

	return reg_morphed_patches;
}