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

	return pattern;
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
	int patch_num = 0;

	for (int i = 0; i < c  ; i++) {
		for (int j = 0; j < r  ; j++) {
			cout << pattern(j,i);
			MatrixXd mask = MatrixXd::Zero(height, width);
			for (int s = columns[i]; s < columns[i + 1]; s++) {
				for (int t = rows[j]; t < rows[j + 1]; t++) {
					mask(t, s) = 255;
				}
			}
			eigen2cv(mask, mask_im);
			imwrite("masks/aligned_mask" + std::to_string(patch_num) + ".png", mask_im);
			patch_num++;
		}
		cout << endl;
	}

	return;
}
