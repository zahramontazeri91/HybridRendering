#include "ImageProcessing.h"
#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <opencv2/core/eigen.hpp>

using Eigen::MatrixXd;
using namespace cv;
using namespace std;

void ControlPoints(Mat movingImage, Mat fixedImage, vector<Point>& movingPoints, vector<Point>& fixedPoints)
{
	int wnd_size = 30;
	int wnd_size_crn = 20;

	Mat edge_moving = EdgeDetector(movingImage);
	Mat edge_fixed = EdgeDetector(fixedImage);
	MatrixXd edge_moving_mat;
	cv2eigen(edge_moving, edge_moving_mat);
	MatrixXd edge_fixed_mat;
	cv2eigen(edge_fixed, edge_fixed_mat);

	//namedWindow("image", CV_WINDOW_AUTOSIZE);
	//imshow("image", fixedImage);

	/// capturing the 4 cornels of aligned mask (fixed image):
	CornerDetector(fixedImage, fixedPoints);
	circle(edge_fixed, fixedPoints[0], 3.0, Scalar(225, 0, 0), 2, 8);
	circle(edge_fixed, fixedPoints[1], 3.0, Scalar(225, 0, 0), 2, 8);
	circle(edge_fixed, fixedPoints[2], 3.0, Scalar(225, 0, 0), 2, 8);
	circle(edge_fixed, fixedPoints[3], 3.0, Scalar(225, 0, 0), 2, 8);

	///Now lets find the corresponding corners in moving image
	//CornerDetector(movingImage, movingPoints);
	// TO DO: lets check Harris if we didn't get one, we can do these:
	MatrixXd movingImg_mat;
	cv2eigen(movingImage, movingImg_mat);

	int min_x_align = fixedPoints[0].x;
	int min_y_align = fixedPoints[0].y;
	int max_x_align = fixedPoints[3].x;
	int max_y_align = fixedPoints[3].y;
	int cnt_moving = 0;

	/// top left corner
	Point crn_tl;
	int min_x = 1000;
	int min_y = 1000;
	for (int i = (min_y_align - wnd_size_crn); i < (min_y_align + wnd_size_crn); i++) {
		for (int j = (min_x_align - wnd_size_crn); j < (min_x_align + wnd_size_crn); j++) {
			if (edge_moving_mat(i, j)) {
				if (i <= min_y && j <= min_x) {
					min_y = i;
					min_x = j;
				}
			}
		}
	}
	crn_tl.x = min_x;
	crn_tl.y = min_y;
	movingPoints.push_back(crn_tl);
	
	/// top right corner
	Point crn_tr;
	int max_x = 1;
	min_y = 1000;
	for (int i = (min_y_align - wnd_size_crn); i < (min_y_align + wnd_size_crn); i++) {
		for (int j = (max_x_align - wnd_size_crn); j < (max_x_align + wnd_size_crn); j++) {
			if (edge_moving_mat(i, j)) {
				if (i <= min_y && j > max_x) {
					min_y = i;
					max_x = j;
				}
			}
		}
	}
	crn_tr.x = max_x;
	crn_tr.y = min_y;
	movingPoints.push_back(crn_tr);

	/// bottom left corner
	Point crn_bl;
	min_x = 1000;
	int max_y = 1;
	for (int i = (max_y_align - wnd_size_crn); i < (max_y_align + wnd_size_crn); i++) {
		for (int j = (min_x_align - wnd_size_crn); j < (min_x_align + wnd_size_crn); j++) {
			if (edge_moving_mat(i, j)) {
				if (i >= max_y && j <= min_x) {
					max_y = i;
					min_x = j;
				}
			}
		}
	}
	crn_bl.x = min_x;
	crn_bl.y = max_y;
	movingPoints.push_back(crn_bl);

	/// bottom right corner
	Point crn_br;
	max_x = 1;
	max_y = 1;
	for (int i = (max_y_align - wnd_size_crn); i < (max_y_align + wnd_size_crn); i++) {
		for (int j = (max_x_align - wnd_size_crn); j < (max_x_align + wnd_size_crn); j++) {
			if (edge_moving_mat(i, j)) {
				if (i >= max_y && j >= max_x) {
					max_y = i;
					max_x = j;
				}
			}
		}
	}
	crn_br.x = max_x;
	crn_br.y = max_y;
	movingPoints.push_back(crn_br);

	circle(edge_moving, crn_tl , 3.0, Scalar(225, 0, 0), 2, 8);
	circle(edge_moving, crn_tr, 3.0, Scalar(225, 0, 0), 2, 8);
	circle(edge_moving, crn_bl, 3.0, Scalar(225, 0, 0), 2, 8);
	circle(edge_moving, crn_br, 3.0, Scalar(225, 0, 0), 2, 8);


	/// ****************Now we want to get the middle points between each 2 corners
	int num_cp = 4 ;
	wnd_size = wnd_size - 1; //it hits negative number for the patches in the first col
	///***********top side of moving image
	int dist = movingPoints[1].x - movingPoints[0].x;
	int step = dist / (num_cp + 1);
	int cnt = 0;
	for (int i = (movingPoints[0].x +1)  ; i < (movingPoints[1].x -1) ; i++) {
		for (int j = (movingPoints[0].y - wnd_size); j < (movingPoints[0].y + wnd_size); j++) {
			if (edge_moving_mat(j,i)) {
				if ((i - movingPoints[0].x +1) % step == 0 && cnt<num_cp ) {
					cnt++;
					movingPoints.push_back(Point(i, j));
					circle(edge_moving, Point(i, j), 3.0, Scalar(255, 0, 0), 2, 8);
					break;
				}
			}
		}
	}
	if (cnt != num_cp) cerr << "Number of moving points is invalid in top side"<<endl;

	/// top side of fixed image	
	dist = fixedPoints[1].x - fixedPoints[0].x;
	step = dist / (num_cp + 1);
	cnt = 0;
	for (int i = (fixedPoints[0].x + 1); i < (fixedPoints[1].x - 1); i++) {
		for (int j = (fixedPoints[0].y - wnd_size); j < (fixedPoints[0].y + wnd_size); j++) {
			if (edge_fixed_mat(j, i)) {
				if ( (i - fixedPoints[0].x + 1) % step == 0 && cnt<num_cp) {
					cnt++;
					fixedPoints.push_back(Point(i, j));
					circle(edge_fixed, Point(i, j), 3.0, Scalar(255, 0, 0), 2, 8);
					break;
				}
			}
		}
	}
	if (cnt != num_cp) cerr << "Number of fixed points is invalid in top side " << endl;

	///***********bottom side of moving image
	dist = movingPoints[3].x - movingPoints[2].x;
	step = dist / (num_cp + 1);
	cnt = 0;
	for (int i = (movingPoints[2].x + 1); i < (movingPoints[3].x - 1); i++) {
		for (int j = (movingPoints[2].y + wnd_size); j > (movingPoints[2].y - wnd_size); j--) {
			if (edge_moving_mat(j, i)) {
				if ((i - movingPoints[2].x + 1) % step == 0 && cnt<num_cp) {
					cnt++;
					movingPoints.push_back(Point(i, j));
					circle(edge_moving, Point(i, j), 3.0, Scalar(255, 0, 0), 2, 8);
					break;
				}
			}
		}
	}
	if (cnt != num_cp) cerr << "Number of moving points is invalid in bottom side" << endl;

	/// bottom side of fixed image	
	dist = fixedPoints[3].x - fixedPoints[2].x;
	step = dist / (num_cp + 1);
	cnt = 0;
	for (int i = (fixedPoints[2].x + 1); i < (fixedPoints[3].x - 1); i++) {
		for (int j = (fixedPoints[2].y + wnd_size); j > (fixedPoints[2].y - wnd_size); j--) {
			if (edge_fixed_mat(j, i)) {
				if ((i - fixedPoints[2].x + 1) % step == 0 && cnt<num_cp) {
					cnt++;
					fixedPoints.push_back(Point(i, j));
					circle(edge_fixed, Point(i, j), 3.0, Scalar(255, 0, 0), 2, 8);
					break;
				}
			}
		}
	}
	if (cnt != num_cp) cerr << "Number of fixed points is invalid in bottom side " << endl;


	///***********left side of moving image
	dist = movingPoints[2].y - movingPoints[0].y;
	step = dist / (num_cp + 1);
	cnt = 0;
	for (int j = (movingPoints[0].y + 1); j < (movingPoints[2].y - 1); j++) {
		for (int i = (movingPoints[0].x - wnd_size); i < (movingPoints[0].x + wnd_size); i++) {
			if (edge_moving_mat(j, i)) {
				if ((j - movingPoints[0].y + 1) % step == 0 && cnt<num_cp) {
					cnt++;
					movingPoints.push_back(Point(i, j));
					circle(edge_moving, Point(i, j), 3.0, Scalar(255, 0, 0), 2, 8);
					break;
				}
			}
		}
	}
	if (cnt != num_cp) cerr << "Number of moving points is invalid in left side" << endl;

	/// left side of fixed image	
	dist = fixedPoints[2].y - fixedPoints[0].y;
	step = dist / (num_cp + 1);
	cnt = 0;
	for (int j = (fixedPoints[0].y + 1); j < (fixedPoints[2].y - 1); j++) {
		for (int i = (fixedPoints[0].x - wnd_size); i < (fixedPoints[0].x + wnd_size); i++) {
			if (edge_fixed_mat(j, i)) {
				if ((j - fixedPoints[0].y + 1) % step == 0 && cnt<num_cp) {
					cnt++;
					fixedPoints.push_back(Point(i, j));
					circle(edge_fixed, Point(i, j), 3.0, Scalar(255, 0, 0), 2, 8);
					break;
				}
			}
		}
	}
	if (cnt != num_cp) cerr << "Number of fixed points is invalid in left side " << endl;

	///***********right side of moving image
	dist = movingPoints[3].y - movingPoints[1].y;
	step = dist / (num_cp + 1);
	cnt = 0;
	for (int j = (movingPoints[1].y + 1); j < (movingPoints[3].y - 1); j++) {
		for (int i = (movingPoints[1].x + wnd_size); i > (movingPoints[1].x - wnd_size); i--) {
			if (edge_moving_mat(j, i)) {
				if ((j - movingPoints[1].y + 1) % step == 0 && cnt<num_cp) {
					cnt++;
					movingPoints.push_back(Point(i, j));
					circle(edge_moving, Point(i, j), 3.0, Scalar(255, 0, 0), 2, 8);
					break;
				}
			}
		}
	}
	if (cnt != num_cp) cerr << "Number of moving points is invalid in right side" << endl;

	/// right side of fixed image	
	dist = fixedPoints[3].y - fixedPoints[1].y;
	step = dist / (num_cp + 1);
	cnt = 0;
	for (int j = (fixedPoints[1].y + 1); j < (fixedPoints[3].y - 1); j++) {
		for (int i = (fixedPoints[1].x + wnd_size); i > (fixedPoints[1].x - wnd_size); i--) {
			if (edge_fixed_mat(j, i)) {
				if ((j - fixedPoints[1].y + 1) % step == 0 && cnt<num_cp) {
					cnt++;
					fixedPoints.push_back(Point(i, j));
					circle(edge_fixed, Point(i, j), 3.0, Scalar(255, 0, 0), 2, 8);
					break;
				}
			}
		}
	}
	if (cnt != num_cp) cerr << "Number of fixed points is invalid in right side " << endl;

	//cout << fixedPoints.size() << endl << movingPoints.size() << endl;

	imshow("edge_moving", edge_moving);
	imshow("edge_fixed", edge_fixed);
	return;
}
