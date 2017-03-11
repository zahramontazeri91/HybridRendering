#ifndef IMAGE_PROCESSING
#define IMAGE_PROCESSING

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <string>

using namespace cv;
using std::string;

Mat ZeroPadding(Mat src, int margin);
Mat EdgeDetector(Mat gray);
void CornerDetector(Mat src_gray, std::vector<Point>& points);
string type2str(int type);
#endif
