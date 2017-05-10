#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/eigen.hpp>
#include "input.h"

using Eigen::MatrixXd;
using namespace cv;
using namespace std;

int patch = 1;

void segmentation() {

	MatrixXd pattern = input_pattern();
	VectorXd columns = input_columns();
	VectorXd rows = input_rows();
	int c = columns.size() - 1;  //6
	int r = rows.size() - 1; //7
	int width = columns[c];
	int height = rows[r];
	Mat mask_im;
	MatrixXd last_half_patch = get_last_half_patch();
	MatrixXd first_half_patch = get_first_half_patch();
	MatrixXd mask = MatrixXd::Zero(height, width);
	int end_col, start_col, end_row, start_row;
	int seg_overlap = 3;

	for (int i = 0; i < c; i++) {
		for (int j = 0; j < r; j++) {
			if (pattern(j, i)) {
				end_col = columns[i + 1];
				start_col = columns[i];
				end_row = rows[j + 1];
				start_row = rows[j];

				if (columns[i + 1] != width)
					end_col = columns[i + 1] + seg_overlap;
				if (columns[i] != 0)
					start_col = columns[i] - seg_overlap;
				if (rows[j + 1] != height)
					end_row = rows[j + 1] + seg_overlap;
				if (rows[j] != 0)
					start_row = rows[j] - seg_overlap;

				for (int s = start_col; s < end_col; s++) {
					for (int t = start_row; t < end_row; t++) {
						mask(t, s) = 255;
					}
				}

				if (last_half_patch(j, i) || rows[j + 1] == height) {
					eigen2cv(mask, mask_im);
					cv::imwrite("input/auto manually mask/patch_" + std::to_string(patch) + ".png", mask_im);
					patch++;
					mask = MatrixXd::Zero(height, width);
				}
			}
		}
	}

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (!pattern(i, j)) {
				end_col = columns[j + 1];
				start_col = columns[j];
				end_row = rows[i + 1];
				start_row = rows[i];

				if (columns[j + 1] != width)
					end_col = columns[j + 1] + seg_overlap;
				if (columns[j] != 0)
					start_col = columns[j] - seg_overlap;
				if (rows[i + 1] != height)
					end_row = rows[i + 1] + seg_overlap;
				if (rows[i] != 0)
					start_row = rows[i] - seg_overlap;

				for (int s = start_col; s < end_col; s++) {
					for (int t = start_row; t < end_row; t++) {
						mask(t, s) = 255;
					}
				}

				if (last_half_patch(i, j) || columns[j + 1] == width) {
					eigen2cv(mask, mask_im);
					cv::imwrite("input/auto manually mask/patch_" + std::to_string(patch) + ".png", mask_im);
					patch++;
					mask = MatrixXd::Zero(height, width);
				}

			}
		}
	}


	return;
}

