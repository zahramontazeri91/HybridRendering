#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/eigen.hpp>

using namespace cv;
using namespace Eigen;
using namespace std; 

int patch_num = 1;

VectorXd input_columns() {

	VectorXd columns(8);
	columns << 0, 50, 134, 227, 290, 360, 440, 457;
	return columns;
}

VectorXd input_rows() {

	VectorXd rows(7);
	rows << 0, 115, 220, 335, 440, 565, 671;
	return rows;
}

MatrixXd input_pattern() {

	MatrixXd pattern(6, 7);
	pattern << ///1 represents warp patch and 0 represents weft patch
		0, 1, 1, 0, 1, 1, 0,
		1, 1, 0, 1, 1, 0, 1,
		1, 0, 1, 1, 0, 1, 1,
		0, 1, 1, 0, 1, 1, 0,
		1, 1, 0, 1, 1, 0, 1,
		1, 0, 1, 1, 0, 1, 1;

	//pattern <<
	//	1, 1, 0, 0, 1, 1, 0,
	//	0, 1, 1, 0, 0, 1, 1,
	//	1, 0, 1, 1, 0, 0, 1,
	//	1, 1, 0, 0, 1, 1, 0,
	//	0, 1, 1, 0, 0, 1, 1,
	//	1, 0, 1, 1, 0, 0, 1;

	return pattern;
}

MatrixXd get_first_half_patch() {

	MatrixXd pattern = input_pattern();
	VectorXd columns = input_columns();
	VectorXd rows = input_rows();
	int c = columns.size() - 1;  //7
	int r = rows.size() - 1; //6
	MatrixXd half_patch_warp(r, c);
	MatrixXd half_patch_weft(r, c);
	MatrixXd first_half_patch(r, c);

	for (int i = 0; i < c; i++) {
		for (int j = 0; j < r; j++) {
			int i_p = i - 1;
			int j_p = j - 1;
			if (j_p < 0) j_p = r - 1;
			if (i_p < 0) i_p = c - 1;

			if (pattern(j, i)) {
				if (pattern(j, i) != pattern(j_p, i))
				{
					first_half_patch(j, i) = 1;
				}
				else
					first_half_patch(j, i) = 0;
			}
			else {
				if (pattern(j,i) != pattern(j, i_p))
				{
					first_half_patch(j,i) = 1;
				}
				else
					first_half_patch(j,i) = 0;
			}
		}
	}

	return first_half_patch;
}
MatrixXd get_last_half_patch() {

	MatrixXd pattern = input_pattern();
	VectorXd columns = input_columns();
	VectorXd rows = input_rows();
	int c = columns.size() - 1;  //7
	int r = rows.size() - 1; //6
	MatrixXd half_patch_warp(r, c);
	MatrixXd half_patch_weft(r, c);
	MatrixXd last_half_patch(r, c);

	for (int i = c-1; i >= 0; i--) {
		for (int j = r-1; j >= 0; j--) {
			int i_p = (i + 1) % c;
			int j_p = (j + 1) % r;

			if (pattern(j, i)) {//warps
				if (pattern(j, i) != pattern(j_p, i))
				{
					last_half_patch(j, i) = 1;
				}
				else
					last_half_patch(j, i) = 0;
			}
			else {//wefts
				if (pattern(j, i) != pattern(j, i_p))
				{
					last_half_patch(j, i) = 1;
				}
				else
					last_half_patch(j, i) = 0;
			}
		}
	}

	return last_half_patch;
}

void get_aligned_masks() {

	MatrixXd pattern = input_pattern();
	VectorXd columns = input_columns();
	VectorXd rows = input_rows();
	int c = columns.size()-1;  //7
	int r = rows.size()-1; //6
	int width = columns[c];
	int height = rows[r];
	Mat mask_im;
	MatrixXd half_patch = get_first_half_patch();
	MatrixXd mask = MatrixXd::Zero(height, width);

	for (int i = 0; i < c  ; i++) {
		for (int j = 0; j < r  ; j++) {
			for (int s = columns[i]; s < columns[i + 1]; s++) {
				for (int t = rows[j]; t < rows[j + 1]; t++) {
					mask(t, s) = 255;
				}
			}
			if (!pattern(j, i) || (pattern(j,i) && (!half_patch(j, i) || rows[j+1]==height) ) ) {
				eigen2cv(mask, mask_im);
				cv::imwrite("input/aligned mask/patch_" + std::to_string(patch_num) + ".png", mask_im);
				patch_num++;
				mask = MatrixXd::Zero(height, width);
			}
		}
	}

	return ;
}

int get_patch_number() {
	return (patch_num - 1);
}
