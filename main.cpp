#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <string>
#include "parser.h"

using namespace cv;
using namespace std;

int main (int argc, char *argv[]) {
	Mat im;
	string file = argv[1];
	im = readfile (480, 640, file);
	cout << "Read file" << endl;
	for (int i = 0 ; i < 4 ; i++) {
		cout << im.at<Vec4f>(240, 320)[i] << " ";
	}
	cout << endl;
}
