#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <string>

#include "arg.h"
#include "parser.h"
#include "method.h"
#include "quilting.h"

using namespace cv;
using namespace std;

int main (int argc, char *argv[]) {
	GeneralParams params;
	if (interpreter (argc, argv, params)) {
		printParameters (params);
		/* Creating a set of images */
		Mat sourceNormalMap, sourceAlbedo, sourceShading;
		Mat targetNormalMap, targetShading;
		Mat newAlbedo, newShading, newImage;
		/* Reading the images specified by the user */
		sourceNormalMap = readfile (480, 640, params.sourceNormalMapFilename);
		sourceAlbedo    = readfile (480, 640, params.sourceAlbedoFilename);
		sourceShading   = readfile (480, 640, params.sourceShadingFilename);
		targetNormalMap = readfile (480, 640, params.targetNormalMapFilename);
		targetShading   = readfile (480, 640, params.targetShadingFilename);
		/* Initializing the other images */
		newAlbedo = Mat(480, 640, CV_32FC(4));
		newShading = Mat(480, 640, CV_32FC(4));
		newImage = Mat(480, 640, CV_32FC(4));
		for (int i = 0 ; i < 480 ; i++) {
			for (int j = 0 ; j < 640 ; j++) {
				for (int k = 0 ; k < 4 ; k++) {
					newAlbedo.at<Vec4f>(i, j)[k] = 0.0f;
					newShading.at<Vec4f>(i, j)[k] = 0.0f;
					newImage.at<Vec4f>(i, j)[k] = 0.0f;
				}
			}
		}
		/* Processing */
		Method *m = new Quilting(sourceNormalMap, sourceAlbedo, sourceShading,
						    	 targetNormalMap, targetShading,
							     newAlbedo, newShading,
								 params.mparams);
		m->textureTransfer();
		m->save();
		delete m;
	}
}

