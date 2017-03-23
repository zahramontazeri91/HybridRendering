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

double percent = 20.0;
int overlap = 5;
int counter = 0;
int width;
int height;
//initialize the theta vector
Eigen::VectorXd theta(4);
double x;
Mat visualization = Mat::zeros(height, width, CV_32FC3);

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


Mat warpRegression(Mat im, int startCol, int endCol, int startRow, int endRow, bool first_half_patch, bool last_half_patch) {

	/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ first half patch
	if (first_half_patch && !last_half_patch) {
		int min_col, max_col, min_row, max_row;
		if (startCol == 0) {
			min_col = startCol;
			min_row = startRow - overlap;
			max_col = endCol + overlap;
			max_row = endRow;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height -1 ) {
				max_row = endRow;
			}
		}
		else if (endCol == width - 1) {
			max_col = endCol;
			min_row = startRow - overlap;
			min_col = startCol - overlap;
			max_row = endRow;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (startRow == 0) {
			min_row = startRow;
			min_col = startCol - overlap;
			max_col = endCol + overlap;
			max_row = endRow;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else if (endRow == height - 1) {
			max_row = endRow;
			min_row = startRow - overlap;
			min_col = startCol - overlap;
			max_col = endCol + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else {
			min_col = startCol - overlap;
			min_row = startRow - overlap;
			max_col = endCol + overlap;
			max_row = endRow;
		}

		MatrixXd im_mat;
		cv2eigen(im, im_mat);

		Point2DVector points;
		int indx, col = 0;

		int transition = int((percent / 100.0) * (endRow - startRow));

		for (int i = startCol; i <= endCol; i++) {
			indx = 0;
			for (int j = startRow; j <= startRow + transition; j++) {
				//circle(visualization, Point(indx, im_mat(j, i)), 1.0, Scalar(255, 0, 0), 2, 8);
				Eigen::Vector2d point;
				point(0) = indx;
				point(1) = im_mat(j, i);
				points.push_back(point);
				indx++;
			}
			col++;
		}

		double d = startRow + transition - min_row;

		MyFunctorNumericalDiff functor;
		functor.Points = points;
		Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

		Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
		//std::cout << "status: " << status << std::endl;

		MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

		//cout << "min_row " << min_row << "startRow " << startRow << "startRow + transition " << startRow + transition << "endRow - transition " << endRow - transition << "endRow " << endRow << "max_row " << max_row << endl;

		for (int i = min_col; i <= max_col; i++) {
			x = 0;
			for (int j = min_row; j < startRow + transition; j++) {
				im_reg_mat(j, i) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
				circle(visualization, Point(x, im_mat(j, i)), 1.0, Scalar(0, 255, 255), 2, 8);
				if (j >= startRow) x++;   // in order to copy the first value for the overlapping region 
			}
		}

		// in order to copy the last value for the constant region in the middle 
		for (int i = min_col; i <= max_col; i++) {
			for (int j = startRow + transition; j <= max_row; j++) {
				im_reg_mat(j, i) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
				circle(visualization, Point(x, im_mat(j, i)), 1.0, Scalar(0, 255, 255), 2, 8);
			}
		}


		Mat regIm, temp;
		cv::eigen2cv(im_reg_mat, temp);
		temp.convertTo(regIm, CV_32FC1);
		//imshow("regularized", regIm);
		//cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
		std::cout << "theta for warp grid_" + std::to_string(counter) << ":  " << endl << theta << std::endl;
		counter++;

		return regIm;
	}

	/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ last half patch
	if (!first_half_patch && last_half_patch) {
		int min_col, max_col, min_row, max_row;
		if (startCol == 0) {
			min_col = startCol;
			min_row = startRow ;
			max_col = endCol + overlap;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (endCol == width - 1) {
			max_col = endCol;
			min_row = startRow;
			min_col = startCol - overlap;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (startRow == 0) {
			min_row = startRow;
			min_col = startCol - overlap;
			max_col = endCol + overlap;
			max_row = endRow + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else if (endRow == height - 1) {
			max_row = endRow;
			min_row = startRow;
			min_col = startCol - overlap;
			max_col = endCol + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else {
			min_col = startCol - overlap;
			min_row = startRow;
			max_col = endCol + overlap;
			max_row = endRow + overlap;
		}

		//cout << min_col << "  " << min_row << "  " << max_col << "  " << max_row << endl;

		MatrixXd im_mat;
		cv2eigen(im, im_mat);
		int transition = int((percent / 100.0) * (endRow - startRow));


		MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

		// in order to copy the last value for the constant region in the middle 
		for (int i = min_col; i <= max_col; i++) {
			for (int j = min_row; j < endRow - transition; j++) {
				im_reg_mat(j, i) = (theta(0)) * tanh( theta(1)* (x + theta(2))) + (theta(3));
				circle(visualization, Point(x, im_mat(j, i)), 1.0, Scalar(0, 255, 255), 2, 8);
			}
		}

		for (int i = min_col; i <= max_col; i++) {
			x = 0;
			for (int j = endRow - transition; j <= max_row; j++) {
				im_reg_mat(j, i) = (theta(0)) * tanh(-1 * theta(1)* (x + theta(2))) + (theta(3));
				circle(visualization, Point(x, im_mat(j, i)), 1.0, Scalar(0, 255, 255), 2, 8);
				if (j <= endRow) x++;
			}
		}

		//cv::imshow("Visualization", visualization);
		cv::imwrite("Reg Visualization/grid_" + std::to_string(counter) + ".png", visualization);
		visualization = Mat::zeros(height, width, CV_32FC3);

		Mat regIm, temp;
		cv::eigen2cv(im_reg_mat, temp);
		temp.convertTo(regIm, CV_32FC1);
		//imshow("regularized", regIm);
		//cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
		std::cout << "theta for warp grid_" + std::to_string(counter) << ":  " << endl << theta << std::endl;
		counter++;

		return regIm;
	}

	/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ neigher patch
	if (!first_half_patch && !last_half_patch) {
		int min_col, max_col, min_row, max_row;
		if (startCol == 0) {
			min_col = startCol;
			min_row = startRow ;
			max_col = endCol + overlap;
			max_row = endRow ;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (endCol == width - 1) {
			max_col = endCol;
			min_row = startRow ;
			min_col = startCol - overlap;
			max_row = endRow ;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (startRow == 0) {
			min_row = startRow;
			min_col = startCol - overlap;
			max_col = endCol + overlap;
			max_row = endRow ;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width ) {
				max_col = endCol;
			}
		}
		else if (endRow == height - 1) {
			max_row = endRow;
			min_row = startRow ;
			min_col = startCol - overlap;
			max_col = endCol + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else {
			min_col = startCol - overlap;
			min_row = startRow ;
			max_col = endCol + overlap;
			max_row = endRow ;
		}


		MatrixXd im_mat;
		cv2eigen(im, im_mat);
		MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

		double y = (theta(0)) * tanh(-1 * theta(1)* (x + theta(2))) + (theta(3));
		// in order to copy the last value for the constant region in the middle 
		for (int i = min_col; i <= max_col; i++) {
			for (int j = min_row; j <= max_row; j++) {
				im_reg_mat(j, i) = y;
			}
		}


		Mat regIm, temp;
		cv::eigen2cv(im_reg_mat, temp);
		temp.convertTo(regIm, CV_32FC1);
		//imshow("regularized", regIm);
		//cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
		std::cout << "theta for warp grid_" + std::to_string(counter) << ":  " << endl << y << std::endl;
		counter++;

		return regIm;
	}

	/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ both patch
	if (first_half_patch && last_half_patch) {
		int min_col, max_col, min_row, max_row;
		if (startCol == 0) {
			min_col = startCol;
			min_row = startRow - overlap;
			max_col = endCol + overlap;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (endCol == width - 1) {
			max_col = endCol;
			min_row = startRow - overlap;
			min_col = startCol - overlap;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (startRow == 0) {
			min_row = startRow;
			min_col = startCol - overlap;
			max_col = endCol + overlap;
			max_row = endRow + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else if (endRow == height - 1) {
			max_row = endRow;
			min_row = startRow - overlap;
			min_col = startCol - overlap;
			max_col = endCol + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
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
		int indx, col = 0;

		int transition = int((percent / 100.0) * (endRow - startRow));

		for (int i = startCol; i <= endCol; i++) {
			indx = 0;
			for (int j = startRow; j <= startRow + transition; j++) {
				//indx = j + col*(endRow - startRow) + 1;
				Eigen::Vector2d point;
				point(0) = indx;
				point(1) = im_mat(j, i);
				points.push_back(point);
				indx++;
			}
			col++;
		}

		double d = startRow + transition - min_row;


		MyFunctorNumericalDiff functor;
		functor.Points = points;
		Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

		Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
		//std::cout << "status: " << status << std::endl;

		MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

		//cout << "min_row " << min_row << "startRow " << startRow << "startRow + transition " << startRow + transition << "endRow - transition " << endRow - transition << "endRow " << endRow << "max_row " << max_row << endl;

		for (int i = min_col; i <= max_col; i++) {
			x = 0;
			for (int j = min_row; j < startRow + transition; j++) {
				im_reg_mat(j, i) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
				if (j >= startRow) x++;   // in order to copy the first value for the overlapping region 
			}
		}

		// in order to copy the last value for the constant region in the middle 
		for (int i = min_col; i <= max_col; i++) {
			for (int j = startRow + transition; j < endRow - transition; j++) {
				im_reg_mat(j, i) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
			}
		}

		for (int i = min_col; i <= max_col; i++) {
			x = 0;
			for (int j = endRow - transition; j <= max_row; j++) {
				im_reg_mat(j, i) = (theta(0)) * tanh(-1 * theta(1)* (x + theta(2))) + (theta(3));
				if (j <= endRow) x++;
			}
		}

		Mat regIm, temp;
		cv::eigen2cv(im_reg_mat, temp);
		temp.convertTo(regIm, CV_32FC1);
		//imshow("regularized", regIm);
		//cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
		std::cout << "theta for warp grid_" + std::to_string(counter) << ":  " << endl << theta << std::endl;
		counter++;

		return regIm;
	}
}

Mat weftRegression(Mat im, int startCol, int endCol, int startRow, int endRow, bool first_half_patch, bool last_half_patch) {

	/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ first half patch
	if (first_half_patch && !last_half_patch) {
		int min_col, max_col, min_row, max_row;
		if (startCol == 0) {
			min_col = startCol;
			min_row = startRow - overlap;
			max_col = endCol;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (endCol == width - 1) {
			max_col = endCol;
			min_row = startRow - overlap;
			min_col = startCol - overlap;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (startRow == 0) {
			min_row = startRow;
			min_col = startCol - overlap;
			max_col = endCol ;
			max_row = endRow + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else if (endRow == height - 1) {
			max_row = endRow;
			min_row = startRow - overlap;
			min_col = startCol - overlap;
			max_col = endCol ;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else {
			min_col = startCol - overlap;
			min_row = startRow - overlap;
			max_col = endCol ;
			max_row = endRow + overlap;
		}
		MatrixXd im_mat;
		cv2eigen(im, im_mat);

		Point2DVector points;
		int indx, row = 0;

		int transition = (percent / 100.0) * (endCol - startCol);
		for (int i = startRow; i <= endRow; i++) {
			indx = 0;
			for (int j = startCol; j <= startCol + transition; j++) {
				//indx = j + row* (endCol - startCol) + 1;
				Eigen::Vector2d point;
				point(0) = indx;
				point(1) = im_mat(i, j);
				points.push_back(point);
				indx++;
			}
			row++;
		}

		MyFunctorNumericalDiff functor;
		functor.Points = points;
		Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

		Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
		//std::cout << "status: " << status << std::endl;

		MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

		double d = startCol + transition - min_col;
		for (int i = min_row; i <= max_row; i++) {
			x = 0;
			for (int j = min_col; j < startCol + transition; j++) {
				//im_reg_mat(i, j) = (theta(0) / 2.0) * tanh( 4.0 / d* (x - d / 2)) + (theta(0) / 2.0);
				im_reg_mat(i, j) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
				if (j >= startCol) x++;	// in order to copy the first value for the overlapping region 
			}
		}

		// in order to copy the last value for the constant region in the middle 
		for (int i = min_row; i <= max_row; i++) {
			for (int j = startCol + transition; j <= max_col; j++) {
				im_reg_mat(i, j) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
			}
		}

		Mat regIm, temp;
		cv::eigen2cv(im_reg_mat, temp);
		temp.convertTo(regIm, CV_32FC1);
		//imshow("regularized", regIm);
		//cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
		std::cout << "theta for weft grid_" + std::to_string(counter) << ":  " << endl << theta << std::endl;
		counter++;

		return regIm;
	}

	/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ last half patch
	if (!first_half_patch && last_half_patch) {
		int min_col, max_col, min_row, max_row;
		if (startCol == 0) {
			min_col = startCol;
			min_row = startRow - overlap;
			max_col = endCol + overlap;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (endCol == width - 1) {
			max_col = endCol;
			min_row = startRow - overlap;
			min_col = startCol ;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (startRow == 0) {
			min_row = startRow ;
			min_col = startCol;
			max_col = endCol + overlap;
			max_row = endRow + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else if (endRow == height - 1) {
			max_row = height;
			min_row = startRow - overlap;
			min_col = startCol ;
			max_col = endCol + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else {
			min_col = startCol ;
			min_row = startRow - overlap;
			max_col = endCol + overlap;
			max_row = endRow + overlap;
		}

		MatrixXd im_mat;
		cv2eigen(im, im_mat);

		Point2DVector points;
		int indx, row = 0;

		int transition = (percent / 100.0) * (endCol - startCol);
		for (int i = startRow; i <= endRow; i++) {
			indx = 0;
			for (int j = startCol; j <= startCol + transition; j++) {
				//indx = j + row* (endCol - startCol) + 1;
				Eigen::Vector2d point;
				point(0) = indx;
				point(1) = im_mat(i, j);
				points.push_back(point);
				indx++;
			}
			row++;
		}


		MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

		double d = startCol + transition - min_col;


		// in order to copy the last value for the constant region in the middle 
		for (int i = min_row; i <= max_row; i++) {
			for (int j = min_col; j < endCol - transition; j++) {
				im_reg_mat(i, j) = (theta(0)) * tanh(-1 * theta(1)* (x + theta(2))) + (theta(3));
			}
		}

		for (int i = min_row; i <= max_row; i++) {
			x = 0;
			for (int j = endCol - transition; j <= max_col; j++) {
				im_reg_mat(i, j) = (theta(0)) * tanh( theta(1)* (x + theta(2))) + (theta(3));
				if (j <= endCol) x++;	// in order to copy the first value for the overlapping region 
			}
		}

		Mat regIm, temp;
		cv::eigen2cv(im_reg_mat, temp);
		temp.convertTo(regIm, CV_32FC1);
		//imshow("regularized", regIm);
		//cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
		std::cout << "theta for weft grid_" + std::to_string(counter) << ":  " << endl << theta << std::endl;
		counter++;

		return regIm;
	}

	/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ neigher patch
	if (!first_half_patch && !last_half_patch) {
		int min_col, max_col, min_row, max_row;
		if (startCol == 0) {
			min_col = startCol;
			min_row = startRow - overlap;
			max_col = endCol ;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (endCol == width - 1) {
			max_col = endCol;
			min_row = startRow - overlap;
			min_col = startCol ;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (startRow == 0) {
			min_row = startRow;
			min_col = startCol ;
			max_col = endCol ;
			max_row = endRow + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else if (endRow == height - 1) {
			max_row = endRow;
			min_row = startRow - overlap;
			min_col = startCol ;
			max_col = endCol ;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else {
			min_col = startCol ;
			min_row = startRow - overlap;
			max_col = endCol;
			max_row = endRow + overlap;
		}

		MatrixXd im_mat;
		cv2eigen(im, im_mat);
		MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

		// in order to copy the last value for the constant region in the middle 
		for (int i = min_row; i <= max_row; i++) {
			for (int j = min_col; j <= max_col; j++) {
				im_reg_mat(i, j) = (theta(0)) * tanh(-1 * theta(1)* (x + theta(2))) + (theta(3));
			}
		}

		Mat regIm, temp;
		cv::eigen2cv(im_reg_mat, temp);
		temp.convertTo(regIm, CV_32FC1);
		//imshow("regularized", regIm);
		//cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
		std::cout << "theta for weft grid_" + std::to_string(counter) << ":  " << endl << 170 << std::endl;
		counter++;

		return regIm;
	}

	/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ both patch
	if (first_half_patch && last_half_patch) {
		int min_col, max_col, min_row, max_row;
		if (startCol == 0) {
			min_col = startCol;
			min_row = startRow - overlap;
			max_col = endCol + overlap;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (endCol == width - 1) {
			max_col = endCol;
			min_row = startRow - overlap;
			min_col = startCol - overlap;
			max_row = endRow + overlap;

			if (startRow == 0) {
				min_row = startRow;
			}
			if (endRow == height - 1) {
				max_row = endRow;
			}
		}
		else if (startRow == 0) {
			min_row = startRow;
			min_col = startCol - overlap;
			max_col = endCol + overlap;
			max_row = endRow + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
				max_col = endCol;
			}
		}
		else if (endRow == height - 1) {
			max_row = endRow;
			min_row = startRow - overlap;
			min_col = startCol - overlap;
			max_col = endCol + overlap;

			if (startCol == 0) {
				min_col = startCol;
			}
			if (endCol == width - 1) {
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

		int transition = (percent / 100.0) * (endCol - startCol);
		for (int i = startRow; i <= endRow; i++) {
			indx = 0;
			for (int j = startCol; j <= startCol + transition; j++) {
				//indx = j + row* (endCol - startCol) + 1;
				Eigen::Vector2d point;
				point(0) = indx;
				point(1) = im_mat(i, j);
				points.push_back(point);
				indx++;
			}
			row++;
		}

		MyFunctorNumericalDiff functor;
		functor.Points = points;
		Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

		Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
		//std::cout << "status: " << status << std::endl;

		MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());

		double d = startCol + transition - min_col;
		for (int i = min_row; i <= max_row; i++) {
			x = 0;
			for (int j = min_col; j < startCol + transition; j++) {
				//im_reg_mat(i, j) = (theta(0) / 2.0) * tanh( 4.0 / d* (x - d / 2)) + (theta(0) / 2.0);
				im_reg_mat(i, j) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
				if (j >= startCol) x++;	// in order to copy the first value for the overlapping region 
			}
		}

		// in order to copy the last value for the constant region in the middle 
		for (int i = min_row; i <= max_row; i++) {
			for (int j = startCol + transition; j < endCol - transition; j++) {
				im_reg_mat(i, j) = (theta(0)) * tanh(theta(1)* (x + theta(2))) + (theta(3));
			}
		}

		for (int i = min_row; i <= max_row; i++) {
			x = 0;
			for (int j = endCol - transition; j <= max_col; j++) {
				im_reg_mat(i, j) = (theta(0)) * tanh(-1 * theta(1)* (x + theta(2))) + (theta(3));
				if (j <= endCol) x++;	// in order to copy the first value for the overlapping region 
			}
		}

		Mat regIm, temp;
		cv::eigen2cv(im_reg_mat, temp);
		temp.convertTo(regIm, CV_32FC1);
		//imshow("regularized", regIm);
		//cv::imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
		std::cout << "theta for weft grid_" + std::to_string(counter) << ":  " << endl << theta << std::endl;
		counter++;

		return regIm;
	}
}

vector<Mat> regularization(vector<Mat> morphed_patches, int padding);
vector<Mat> regularization(int padding ) {
	vector<Mat> morphed_patches;
	for (int i = 0; i < 30; i++) {
		morphed_patches.push_back(imread("Morphed Patches/morphed_patch_" + std::to_string(i) + ".png", CV_LOAD_IMAGE_GRAYSCALE));
	}
	vector<Mat> reg_morphed_patches = regularization(morphed_patches, padding);
	return reg_morphed_patches;
}
vector<Mat> regularization(vector<Mat> morphed_patches, int padding) {

	vector<Mat> reg_morphed_patches;
	int patch_num = get_patch_number();
	MatrixXd pattern = input_pattern();
	MatrixXd first_half_patch = get_first_half_patch();
	MatrixXd last_half_patch = get_last_half_patch();
	VectorXd columns = input_columns();
	VectorXd rows = input_rows();
	int cs = pattern.cols();  //7
	int rs = pattern.rows(); //6
	width = columns[cs];
	height = rows[rs];
	theta << 170.0 / 2.0, 4.0 / 44, -1.0 * (44 / 2.0), 170.0 / 2.0;

	//int startCol, int endCol, int startRow, int endRow,
	int i = 0;
	Mat temp = Mat(height, width, CV_32FC1, cvScalar(0.));
	for (int c = 0; c < cs; c++) {
		for (int r = 0; r < rs; r++) {
			if (pattern(r, c)) {
				Mat temp2 = warpRegression(morphed_patches[i], columns[c], columns[c + 1] -1 , rows[r], rows[r + 1] -1 , first_half_patch(r, c), last_half_patch(r, c));
				temp = temp + temp2;
				if (last_half_patch(r, c) || rows[r + 1] == height) {
					reg_morphed_patches.push_back(temp);
					i++;
					temp = Mat(height, width, CV_32FC1, cvScalar(0.));
				}
			}
			else if (!pattern(r, c) ) {
				Mat temp2 = weftRegression(morphed_patches[i], columns[c], columns[c + 1] - 1, rows[r], rows[r + 1] - 1, first_half_patch(r, c), last_half_patch(r, c));
				temp = temp + temp2;
				if (last_half_patch(r, c) || columns[c + 1] == width) {
					reg_morphed_patches.push_back(temp);
					i++;
					temp = Mat(height, width, CV_32FC1, cvScalar(0.));
				}
			}
			/// go to next grid if it is the last_half_patch

		}
	}

	return reg_morphed_patches;
}