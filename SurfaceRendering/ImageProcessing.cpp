#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <vector>

using namespace cv;

/// Global Variables
Mat  dst;
Scalar value;
//char* window_name = "copyMakeBorder Demo";
RNG rng(12345);

/// Global variables for CornerDetector

int thresh = 200;
int max_thresh = 255;
char* source_window = "Source image";
char* corners_window = "Corners detected";

Mat ZeroPadding(Mat src)
{
	dst = src;
	value = Scalar(0, 0, 0);
	copyMakeBorder(src, dst, 30, 30, 30, 30, BORDER_CONSTANT, value);

	// Create window
	//namedWindow(window_name, CV_WINDOW_AUTOSIZE);
	//imshow(window_name, dst);

	return dst;
}

Mat EdgeDetector(Mat gray)
{
	Mat edge, draw;
	Canny(gray, edge, 50, 150, 3);
	edge.convertTo(draw, CV_8U);
	return draw;
}

void CornerDetector(Mat src_gray, std::vector<Point>& points)
{

	Mat dst, dst_norm, dst_norm_scaled;
	dst = Mat::zeros(src_gray.size(), CV_32FC1);

	/// Detector parameters
	int blockSize = 2;
	int apertureSize = 3;
	double k = 0.04;

	/// Detecting corners
	cornerHarris(src_gray, dst, blockSize, apertureSize, k, BORDER_DEFAULT);

	/// Normalizing
	normalize(dst, dst_norm, 0, 255, NORM_MINMAX, CV_32FC1, Mat());
	convertScaleAbs(dst_norm, dst_norm_scaled);

	/// Drawing a circle around corners
	for (int j = 0; j < dst_norm.rows; j++)
	{
		for (int i = 0; i < dst_norm.cols; i++)
		{
			if ((int)dst_norm.at<float>(j, i) > thresh)
			{
				circle(dst_norm_scaled, Point(i, j), 5, Scalar(0), 2, 8, 0);
				points.push_back(Point(i, j));
			}
		}
	}
	/// Showing the result
	namedWindow(corners_window, CV_WINDOW_AUTOSIZE);
	imshow(corners_window, dst_norm_scaled);
}