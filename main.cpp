#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <string>

#include "parser.h"
#include "arg.h"

using namespace cv;
using namespace std;

int main (int argc, char *argv[]) {
	ReconstructionParams params;
	if (interpreter (argc, argv, params)) {
		Mat im;
		string file = params.sourceName;
		file += "_n.txt";
		im = readfile (480, 640, file);
		cout << "Read file" << endl;
		for (int i = 0 ; i < 4 ; i++) {
			cout << im.at<Vec4f>(240, 320)[i] << " ";
		}
		cout << endl;
	}
	printParameters (params);
}
