#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "method.h"

#define clamp(a, x, b) ((x) <= (a) ? (a) : ((b) <= (x) ? (b) : (x)))


Method::Method (Mat & _sourceNormalMap,
				Mat & _sourceAlbedo,
				Mat & _sourceShading,
				Mat & _targetNormalMap,
				Mat & _targetShading,
				Mat & _newAlbedo,
				Mat & _newShading,
				MethodParams & _params) :
	sourceNormalMap (_sourceNormalMap),
	sourceAlbedo (_sourceAlbedo),
	sourceShading (_sourceShading),
	targetNormalMap (_targetNormalMap),
	targetShading (_targetShading),
	newAlbedo (_newAlbedo),
	newShading (_newShading),
	params (_params)
{

}


Method::~Method () {

}


void Method::save () {
	_save("./sourceNormalMap.png", sourceNormalMap);
	_save("./sourceAlbedo.png", sourceAlbedo);
	_save("./targetNormalMap.png", targetNormalMap);
	_save("./target.png", newAlbedo);
}


void Method::_save (string filename, Mat m) {
	Mat x (480, 640, CV_8UC4);
	for (int i = 0 ; i < 480 ; i++) {
		for (int j = 0 ; j < 640 ; j++) {
			Vec4b & bgra = x.at<Vec4b>(i, j);
			bgra[0] = int(clamp(0, 255 * m.at<Vec4f>(i, j)[2], 255));
			bgra[1] = int(clamp(0, 255 * m.at<Vec4f>(i, j)[1], 255));
			bgra[2] = int(clamp(0, 255 * m.at<Vec4f>(i, j)[0], 255));
			bgra[3] = int(clamp(0, 255 * m.at<Vec4f>(i, j)[3], 255));
		}
	}
	imwrite(filename, x);
}
