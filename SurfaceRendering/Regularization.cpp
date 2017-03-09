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

int counter = 0;
//Point2DVector GeneratePoints();

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
			fvec(i) = this->Points[i](1) - (x(0) * sin(x(1)*this->Points[i](0) + x(2)) + x(3));
		}

		return 0;
	}

	Point2DVector Points;

	int inputs() const { return 4; } // There are two parameters of the model
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


Mat warpRegression(Mat im, int startCol, int endCol) {

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
	Eigen::VectorXd theta(4);
	theta << 233/2, 2*3.14/335, 0, 0;
	//x.fill(4.0f);

	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
	std::cout << "status: " << status << std::endl;

	MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());
	for (int i = startCol; i <= endCol; i++) {
		for (int j = 0; j < im_mat.rows(); j++) {
			double x = j;
			double y = theta(0)*sin(theta(1)*x + theta(2)) + theta(3);
			cout << y << endl;
			im_reg_mat(j, i) = y;
		}
	}


	Mat regIm, temp;
	eigen2cv(im_reg_mat, temp);
	temp.convertTo(regIm, CV_8UC1);
	imshow("regularized", regIm);
	imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
	counter++;

	std::cout << "x that minimizes the function: " << std::endl << theta << std::endl;
	return regIm;
}

Mat weftRegression(Mat im, int startRow, int endRow) {

	MatrixXd im_mat;
	cv2eigen(im, im_mat);

	Point2DVector points;
	int indx, row = 0;
	//TO DO: make sure you changed the size of image back
	for (int i = startRow; i <= endRow; i++) {
		for (int j = 0; j < im_mat.cols(); j++) {
			indx = j + row*im_mat.cols();
			Eigen::Vector2d point;
			point(0) = indx;
			point(1) = im_mat(i,j);
			points.push_back(point);

		}
		row++;
	}

	//unsigned int numberOfPoints = 200;
	//Point2DVector points = GeneratePoints(numberOfPoints);

	//initialize the theta vector
	Eigen::VectorXd theta(4);
	theta << 233 / 2, 2 * 3.14 / 335, 0, 0;
	//x.fill(4.0f);

	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);

	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(theta);
	std::cout << "status: " << status << std::endl;

	MatrixXd im_reg_mat = MatrixXd::Zero(im_mat.rows(), im_mat.cols());
	for (int i = startRow; i <= endRow; i++) {
		for (int j = 0; j < im_mat.cols(); j++) {
			double x = j;
			double y = theta(0)*sin(theta(1)*x + theta(2)) + theta(3);
			cout << y << endl;
			im_reg_mat(i,j) = y;
		}
	}


	Mat regIm, temp;
	eigen2cv(im_reg_mat, temp);
	temp.convertTo(regIm, CV_8UC1);
	imshow("regularized", regIm);
	imwrite("Regressioned Patches/reg_patch_" + std::to_string(counter) + ".png", regIm);
	counter++;

	std::cout << "x that minimizes the function: " << std::endl << theta << std::endl;
	return regIm;
}

void regularization(vector<Mat> morphed_patches, vector<Mat>& warpsReg, vector<Mat>& weftsReg);
void regularization(vector<Mat>& warpsReg, vector<Mat>& weftsReg) {
	vector<Mat> morphed_patches;
	for (int i = 0; i < 30; i++) {
		morphed_patches.push_back(imread("Morphed Patches/morphed_patch_" + std::to_string(i) + ".png", CV_LOAD_IMAGE_GRAYSCALE));
	}
	regularization(morphed_patches, warpsReg, weftsReg);
	return;
}
void regularization(vector<Mat> morphed_patches, vector<Mat>& warpsReg, vector<Mat>& weftsReg) {

	///putting together the patches from same yarn
	vector<Mat> warps;
	vector<Mat> wefts;

	warps.push_back(morphed_patches[1] + morphed_patches[3]);
	warps.push_back(morphed_patches[4] + morphed_patches[6]);
	warps.push_back(morphed_patches[8] + morphed_patches[10] + morphed_patches[12]);
	warps.push_back(morphed_patches[14] + morphed_patches[16]);
	warps.push_back(morphed_patches[17] + morphed_patches[19]);
	warps.push_back(morphed_patches[21] + morphed_patches[23] + morphed_patches[25]);
	warps.push_back(morphed_patches[27] + morphed_patches[29]);

	wefts.push_back(morphed_patches[0] + morphed_patches[13] + morphed_patches[26]);
	wefts.push_back(morphed_patches[9] + morphed_patches[22]);
	wefts.push_back(morphed_patches[5] + morphed_patches[18]);
	wefts.push_back(morphed_patches[2] + morphed_patches[15] + morphed_patches[28]);
	wefts.push_back(morphed_patches[11] + morphed_patches[24]);
	wefts.push_back(morphed_patches[7] + morphed_patches[20]);

	///Apply the regression of each yarn
	weftsReg.push_back(weftRegression(wefts[0], 0, 115));
	weftsReg.push_back(weftRegression(wefts[1], 115, 220));
	weftsReg.push_back(weftRegression(wefts[2], 220, 335));
	weftsReg.push_back(weftRegression(wefts[3], 335, 440));
	weftsReg.push_back(weftRegression(wefts[4], 440, 565));
	weftsReg.push_back(weftRegression(wefts[5], 565, 670));

	warpsReg.push_back(warpRegression(warps[0], 0 , 50));
	warpsReg.push_back(warpRegression(warps[1], 50, 134));
	warpsReg.push_back(warpRegression(warps[2], 134,227));
	warpsReg.push_back(warpRegression(warps[3], 227, 290));
	warpsReg.push_back(warpRegression(warps[4], 290, 360));
	warpsReg.push_back(warpRegression(warps[5], 360, 440));
	warpsReg.push_back(warpRegression(warps[6], 440,456));

	return;
}