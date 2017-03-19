#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/core/eigen.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "ImageProcessing.h"
#include "Regularization.h"
#include "ControlPoints.h"
#include "Morphing.h"
#include "input.h"
#include <time.h>

using Eigen::MatrixXd;
using namespace cv;
using namespace std;

int main(int argc, char** argv)
{
	clock_t init, final;
	init = clock();

	Mat grayImage;
	vector<Mat> alignedMasks;
	vector<Mat> alignedMasks_pad;
	vector<Mat> masks;
	vector<Mat> masks_pad;
	vector<Mat> patches;
	vector<Mat> patches_pad;
	vector<Mat> morphed_patches;
	vector<Mat> reg_patches;
	vector<Mat> reg_yarns;
	vector<Mat> unmorphed_reg_patches;
	int padding = 30;
	vector < vector<Point> > movingPoints;
	vector < vector<Point> > fixedPoints;
	vector<Mat> gauss_masks;
	
	/*
	* Input Reading section
	*/ 
	std::cout << "***************************Reading input section **************" << endl;
	get_aligned_masks();
	int patch_num = get_patch_number();

	grayImage = imread("input/height.exr", IMREAD_GRAYSCALE); // Read the height map
	if (!grayImage.data) // Check for invalid input
	{
		std::cout << "Could not open or find the image" << std::endl;
		return -1;
	}


	/*
	* Segmentation section
	*/
	std::cout << "***************************segmentation section **************" << endl;
	for (int i = 0; i < patch_num; i++) {
		alignedMasks.push_back( imread("input/aligned mask/patch_"+ to_string(i+1)+ ".png", CV_LOAD_IMAGE_GRAYSCALE) );
	}

	for (int i = 0; i < patch_num; i++) {
		masks.push_back(imread("input/manually mask/patch_" + to_string(i + 1) + ".png", CV_LOAD_IMAGE_GRAYSCALE));
	}
	///mask the input to get all the patches 
	for (int i = 0; i < patch_num; i++) {
		Mat temp;
		grayImage.copyTo(temp, masks[i]);
		patches.push_back(temp);
		morphed_patches.push_back(Mat::zeros(grayImage.rows, grayImage.cols, CV_32FC1));
		unmorphed_reg_patches.push_back(Mat::zeros(grayImage.rows, grayImage.cols, CV_32FC1));
		imwrite("Patches/patch_" + std::to_string(i) + ".png", temp);
	}


	/*
	* Control Point section
	*/
	cout << "************************ Control Point section ********************* " << endl;
	for (int i = 0; i < patch_num; i++) {
		masks_pad.push_back( ZeroPadding(masks[i], padding));
		alignedMasks_pad.push_back(ZeroPadding(alignedMasks[i], padding));
		patches_pad.push_back(ZeroPadding(patches[i], padding));

		vector<Point> movingPoints_tmp;
		vector<Point> fixedPoints_tmp;
		ControlPoints(masks_pad[i], alignedMasks_pad[i], movingPoints_tmp, fixedPoints_tmp);
		movingPoints.push_back(movingPoints_tmp);
		fixedPoints.push_back(fixedPoints_tmp);
	}


	/*
	* morphing section using the control points for each patch
	*/
	//cout << "************************ morphing section ********************* " << endl;
	//for (int i = 0; i < patch_num; i++) {
	//	Mat temp;
	//	Mat grayImage_pad = ZeroPadding(grayImage, padding);
	//	morphing(temp, grayImage_pad, movingPoints[i], fixedPoints[i]);

	//	///change the size of the results back, before padding
	//	// Setup a rectangle to define the region of interest with width and height
	//	cv::Rect myROI(padding, padding, grayImage.cols, grayImage.rows);
	//	//// Crop the full image to that image contained by the rectangle myROI
	//	cv::Mat croppedTemp = temp(myROI);
	//	croppedTemp.copyTo(morphed_patches[i], alignedMasks[i]);
	//	imwrite("Morphed Patches/morphed_patch_" + std::to_string(i) + ".png", morphed_patches[i]);
	//}


	/*
	* regularization section for each patch
	*/
	std::cout << "************************Regularization section ************** " << endl;
	vector<Mat> reg_morphed_patches;
	reg_morphed_patches = regularization(fixedPoints, padding);


	/*
	* Undo_Morphing section for each patch
	*/
	std::cout << "************************Undo_Morphing section *************" << endl;
	cv::Rect myROI(padding, padding, grayImage.cols, grayImage.rows);

	for (int i = 0; i < patch_num; i++) {
		Mat morphed_yarn;
		Mat reg_patch_pad = ZeroPadding(reg_morphed_patches[i], padding);
		morphing(morphed_yarn, reg_patch_pad, fixedPoints[i], movingPoints[i]);
		///change the size of the results back to before padding
		/// Crop the full image to that image contained by the rectangle myROI
		cv::Mat cropped_morphed_yarn = morphed_yarn(myROI);
		Mat temp;
		GaussianBlur(masks[i], temp, Size(7,7), 0 , 0 );//applying Gaussian filter 
		gauss_masks.push_back(temp);

		///make the type of Mats same and ready for element-wise maltiplication:
		gauss_masks[i].convertTo(gauss_masks[i], CV_32FC1, 1.0 / 255.0);
		cropped_morphed_yarn.convertTo(cropped_morphed_yarn, CV_32FC1, 1.0 / 255.0);

		//string ty = type2str(unmorphed_reg_patches[i].type());
		//printf("Matrix: %s %dx%d \n", ty.c_str(), unmorphed_reg_patches[i].cols, unmorphed_reg_patches[i].rows);
		//string ty2 = type2str(cropped_morphed_yarn.type());
		//printf("Matrix: %s %dx%d \n", ty2.c_str(), cropped_morphed_yarn.cols, cropped_morphed_yarn.rows);
		//string ty3 = type2str(gauss_masks[i].type());
		//printf("Matrix: %s %dx%d \n", ty3.c_str(), gauss_masks[i].cols, gauss_masks[i].rows);

		unmorphed_reg_patches[i] = cropped_morphed_yarn.mul(gauss_masks[i]);
		imwrite("Unmorphed Reg/unmorphed_reg_patch_" + std::to_string(i) + ".png", unmorphed_reg_patches[i]);

	}


	/*
	* Putting the patches together and blending
	*/
	std::cout << "************************Blending section *************" << endl;
	/// TO DO: why it doesn't work if I read from file instead
	//for (int i = 0; i < patch_num; i++) { 
	//	unmorphed_reg_patches[i] = cv::imread("Unmorphed Reg/unmorphed_reg_patch_" + std::to_string(i) + ".png", CV_LOAD_IMAGE_GRAYSCALE);
	//}
	Mat regularized = unmorphed_reg_patches[0];
	for (int i = 1; i < patch_num; i++) {
		regularized += unmorphed_reg_patches[i];
	}
	cv::imshow("Regularized Map", regularized);
	//imwrite("Output/regularized.jpg", regularized);

	/*
	* Obtaining Residual map
	*/
	std::cout << "************************Residual map section *************" << endl;
	Mat residual = Mat::zeros(grayImage.rows, grayImage.cols, CV_32FC1);

	grayImage.convertTo(grayImage, CV_32FC1, 1.0 / 255.0);
	regularized.convertTo(regularized, CV_32FC1, 1.0 / 255.0);

	///make the type of Mats same and ready for element-wise subtraction:
	//string ty3 = type2str(regularized.type());
	//printf("Matrix: %s %dx%d \n", ty3.c_str(), regularized.cols, regularized.rows);
	//string ty4 = type2str(residual.type());
	//printf("Matrix: %s %dx%d \n", ty4.c_str(), residual.cols, residual.rows);
	//string ty5 = type2str(grayImage.type());
	//printf("Matrix: %s %dx%d \n", ty5.c_str(), grayImage.cols, grayImage.rows);

	subtract(grayImage, regularized, residual);
	//residual = grayImage - regularized;

	///get the absolute value:
	double min;
	double max;
	cv::minMaxIdx(residual, &min, &max);
	cv::Mat adjMap;
	// expand your range to 0..255. Similar to histEq();
	residual.convertTo(adjMap, CV_8UC1, 255 / (max - min), -255 * min / (max - min));
	cv::Mat falseColorsMap;
	cv::applyColorMap(adjMap, falseColorsMap, cv::COLORMAP_AUTUMN);
	cv::imshow("Residual Map", residual);

	cv::imshow("Height Map", grayImage);
	//////******************************* Masking section
	///// now let's separate out the patches from regularized yarns
	//Mat temp;
	//weftsReg[0].copyTo(temp, alignedMasks[0]);
	//reg_patches.push_back(temp);
	//warpsReg[0].copyTo(temp, alignedMasks[1]);
	//reg_patches.push_back(temp);
	//weftsReg[3].copyTo(temp, alignedMasks[2]);
	//reg_patches.push_back(temp);
	//warpsReg[0].copyTo(temp, alignedMasks[3]);
	//reg_patches.push_back(temp);
	//warpsReg[1].copyTo(temp, alignedMasks[4]);
	//reg_patches.push_back(temp);
	//weftsReg[2].copyTo(temp, alignedMasks[5]);
	//reg_patches.push_back(temp);
	//warpsReg[1].copyTo(temp, alignedMasks[6]);
	//reg_patches.push_back(temp);
	//weftsReg[5].copyTo(temp, alignedMasks[7]);
	//reg_patches.push_back(temp);
	//warpsReg[2].copyTo(temp, alignedMasks[8]);
	//reg_patches.push_back(temp);
	//weftsReg[1].copyTo(temp, alignedMasks[9]);
	//reg_patches.push_back(temp);
	//warpsReg[2].copyTo(temp, alignedMasks[10]);
	//reg_patches.push_back(temp);
	//weftsReg[4].copyTo(temp, alignedMasks[11]);
	//reg_patches.push_back(temp);
	//warpsReg[2].copyTo(temp, alignedMasks[12]);
	//reg_patches.push_back(temp);
	//weftsReg[0].copyTo(temp, alignedMasks[13]);
	//reg_patches.push_back(temp);
	//warpsReg[3].copyTo(temp, alignedMasks[14]);
	//reg_patches.push_back(temp);
	//weftsReg[4].copyTo(temp, alignedMasks[15]);
	//reg_patches.push_back(temp);
	//warpsReg[3].copyTo(temp, alignedMasks[16]);
	//reg_patches.push_back(temp);
	//warpsReg[4].copyTo(temp, alignedMasks[17]);
	//reg_patches.push_back(temp);
	//weftsReg[2].copyTo(temp, alignedMasks[18]);
	//reg_patches.push_back(temp);
	//warpsReg[4].copyTo(temp, alignedMasks[19]);
	//reg_patches.push_back(temp);
	//weftsReg[5].copyTo(temp, alignedMasks[20]);
	//reg_patches.push_back(temp);
	//warpsReg[5].copyTo(temp, alignedMasks[21]);
	//reg_patches.push_back(temp);
	//weftsReg[1].copyTo(temp, alignedMasks[22]);
	//reg_patches.push_back(temp);
	//warpsReg[5].copyTo(temp, alignedMasks[23]);
	//reg_patches.push_back(temp);
	//weftsReg[4].copyTo(temp, alignedMasks[24]);
	//reg_patches.push_back(temp);
	//warpsReg[5].copyTo(temp, alignedMasks[25]);
	//reg_patches.push_back(temp);
	//weftsReg[0].copyTo(temp, alignedMasks[26]);
	//reg_patches.push_back(temp);
	//warpsReg[6].copyTo(temp, alignedMasks[27]);
	//reg_patches.push_back(temp);
	//weftsReg[3].copyTo(temp, alignedMasks[28]);
	//reg_patches.push_back(temp);
	//warpsReg[6].copyTo(temp, alignedMasks[29]);
	//reg_patches.push_back(temp);

	/// To check FastLaplace() is working
	//MatrixXd Dx(4,5);
	//Dx << 8.9, 8.9, 8.9, 8.9, 8.9,
	//	8.4, 0, 0, 0, 9.2,
	//	7.2, 0, 10.2, 0, 9.4,
	//	6.1, 6.8, 7.7, 8.7, 6.1;
	//MatrixXd Dy(4,5);
	//Dy << 8.9, 8.9, 8.9, 8.9, 8.9,
	//	8.4, 0, 0, 0, 9.2,
	//	7.2, 0, 0, 0, 9.4,
	//	6.1, 6.8, 7.7, 8.7, 6.1;
	//MatrixXd temp(4,5);
	//temp << 8.9, 8.9, 8.9, 8.9, 8.9,
	//	8.4, 0, 0, 0, 9.2,
	//	7.2, 0, 1, 0, 9.4,
	//	6.1, 6.8, 7.7, 8.7, 6.1;
	//MatrixXd Lx(4,5), Ly(4,5);

	final = clock() - init;
	std::cout << "the program is finished in "<< (double)final / ((double)CLOCKS_PER_SEC);

	waitKey(0); // Wait for a keystroke in the window
	return 0;
}

