#ifndef IMAGE_PROCESSING
#define IMAGE_PROCESSING

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;

Mat ZeroPadding(Mat src, int margin);
Mat EdgeDetector(Mat gray);
void CornerDetector(Mat src_gray, std::vector<Point>& points);
#endif
