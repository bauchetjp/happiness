#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <string>
#include <cstdio>

using namespace cv;
using namespace std;

Mat readfile (int height, int width, string name) {
	Mat im (height, width, CV_32FC(4));
	FILE* imStream = fopen(name.c_str(), "r");
	if (imStream == NULL) {
		cout << "Warning : failed to open file " << name << endl;
		for (int i = 0 ; i < im.rows ; i++) {
			for (int j = 0 ; j < im.cols ; j++) {
				im.at<Vec4f>(i, j)[0] = 0.0f;
				im.at<Vec4f>(i, j)[1] = 0.0f;
				im.at<Vec4f>(i, j)[2] = 0.0f;
				im.at<Vec4f>(i, j)[3] = 0.0f;
			}
		}
	} else {
		float r, g, b, a;
		for (int i = 0 ; i < im.rows ; i++) {
			for (int j = 0 ; j < im.cols ; j++) {
				fscanf(imStream, "%f %f %f %f\n", &r, &g, &b, &a);
				im.at<Vec4f>(i, j)[0] = r;
				im.at<Vec4f>(i, j)[1] = g;
				im.at<Vec4f>(i, j)[2] = b;
				im.at<Vec4f>(i, j)[3] = a;
			}
		}
		fclose(imStream);		
	}
	return im;
}

