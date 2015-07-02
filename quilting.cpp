#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/flann/flann.hpp>
#include <list>
#include <vector>
#include <utility>
#include <iostream>
#include "method.h"
#include "quilting.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define clamp(a, x, b) ((x) <= (a) ? (a) : ((b) <= (x) ? (b) : (x)))

using namespace cv;
using namespace std;

Quilting::Quilting (Mat & _sourceNormalMap,
		   			Mat & _sourceAlbedo,
	 			    Mat & _sourceShading,
		   			Mat & _targetNormalMap,
		   			Mat & _targetShading,
		   			Mat & _newAlbedo,
		   			Mat & _newShading,
					MethodParams & _params) :
	Method (_sourceNormalMap,
			_sourceAlbedo,
			_sourceShading,
			_targetNormalMap,
			_targetShading,
			_newAlbedo,
			_newShading,
			_params)
{
	centers = list<Point2i>();
}


Quilting::~Quilting() {
	centers.clear();
	features.release();
}


void Quilting::textureTransfer () {
	for (int i = 0 ; i < params.quiltingIterations ; i++) {
		if (i == 0) {
			radius = params.quiltingRadius;
			overlap = max(2, radius / 3);	
			dimensionality = 3 * ((2 * radius + 1) * (2 * radius + 1) + 2 * overlap * (2 * radius + 1) - overlap * overlap);
		} else {
			radius = max(2, (2 * radius) / 3);
			overlap = 0;
			dimensionality = 3 * (8 * radius * radius + 8 * radius + 2);
		}
		centers.clear();
		features.release();
		entries = detectFeatures ();
		cout << "Iteration " << i << ", features : " << entries << ", dimensionality : " << dimensionality << endl;
		if (params.quiltingIterations == 1) {
			alpha = 0.9;
		} else {
			alpha = 0.9 - 0.8 * float(i) / (params.quiltingIterations - 1);
		}
		beta = 1 - alpha;
		buildFeatures (i);
		matchFeatures (i);
	}
}


int Quilting::detectFeatures () {
	for (int i = 0 ; i < sourceNormalMap.rows ; i++) {
		for (int j = 0 ; j < sourceNormalMap.cols ; j++) {
			bool fullyTransparent = true;
			for (int y = max(0, i - radius) ; y <= min(sourceNormalMap.rows - 1, i + radius) && fullyTransparent ; y++) {
				for (int x = max(0, j - radius) ; x <= min(sourceNormalMap.cols - 1, j + radius) && fullyTransparent ; x++) {
					if (sourceNormalMap.at<Vec4f>(y, x)[3] > 0.9) {
						centers.push_back(Point2i(j, i));
						fullyTransparent = false;
					}
				}
			}
		}
	}
	return centers.size();
}


void Quilting::buildFeatures (int currentIteration) {
	features.create(entries, dimensionality, CV_32F);
	list<Point2i>::iterator it;
	int address = 0;
	int shift = 3 * (2 * radius + 1) * (2 * radius + 1);
	int shift_2 = 3 * overlap * (2 * radius + 1);
	for (it = centers.begin() ; it != centers.end() ; it++) {
		int i = it->y;
		int j = it->x;
		/* Copying the normal map in the feature */
		for (int y = i - radius ; y <= i + radius ; y++) {
			for (int x = j - radius ; x <= j + radius ; x++) {
				int yp = y - (i - radius);
				int xp = x - (j - radius);
				/* Checks if the current point is inside the normal map */
				if (y >= 0 && y < sourceNormalMap.rows && x >= 0 && x < sourceNormalMap.cols) {
					features.at<float>(address, 3 * (yp * ((radius << 1) + 1) + xp) + 0) = sqrt(alpha) * sourceNormalMap.at<Vec4f>(y, x)[0];
					features.at<float>(address, 3 * (yp * ((radius << 1) + 1) + xp) + 1) = sqrt(alpha) * sourceNormalMap.at<Vec4f>(y, x)[1];
					features.at<float>(address, 3 * (yp * ((radius << 1) + 1) + xp) + 2) = sqrt(alpha) * sourceNormalMap.at<Vec4f>(y, x)[2];
				} else {
					features.at<float>(address, 3 * (yp * ((radius << 1) + 1) + xp) + 0) = 0.0f;
					features.at<float>(address, 3 * (yp * ((radius << 1) + 1) + xp) + 1) = 0.0f;
					features.at<float>(address, 3 * (yp * ((radius << 1) + 1) + xp) + 2) = 0.0f;
				}
			}
		}
		/* Copying the albedo in the feature */
		if (currentIteration == 0) {
			/* We copy the top left corner of the tile */
			for (int y = i - radius ; y <= i + radius ; y++) {
				for (int x = j - radius ; x <= j + radius ; x++) {
					int yp = y - (i - radius);
					int xp = x - (j - radius);
					/* Checks if the current point is inside the normal map */
					if (y >= 0 && y < sourceNormalMap.rows && x >= 0 && x < sourceNormalMap.cols) {
						if (yp < overlap) {
							features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 0) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[0];
							features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 1) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[1];
							features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 2) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[2];
						} else if (xp < overlap) {
							features.at<float>(address, shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 0) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[0];
							features.at<float>(address, shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 1) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[1];
							features.at<float>(address, shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 2) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[2];
						}
					} else {
						if (yp < overlap) {
							features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 0) = 0.0f;
							features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 1) = 0.0f;
							features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 2) = 0.0f;
						} else if (xp < overlap) {
							features.at<float>(address, shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 0) = 0.0f;
							features.at<float>(address, shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 1) = 0.0f;
							features.at<float>(address, shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 2) = 0.0f;
						}
					}
				}
			}
		} else {
			/* We copy the entire tile */
			for (int y = i - radius ; y <= i + radius ; y++) {
				for (int x = j - radius ; x <= j + radius ; x++) {
					int yp = y - (i - radius);
					int xp = x - (j - radius);
					if (y >= 0 && y < sourceNormalMap.rows && x >= 0 && x < sourceNormalMap.cols) {
						features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 0) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[0];
						features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 1) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[1];
						features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 2) = sqrt(beta) * sourceAlbedo.at<Vec4f>(y, x)[2];
					} else {
						features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 0) = 0.0f;
						features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 1) = 0.0f;
						features.at<float>(address, shift + 3 * (yp * ((radius << 1) + 1) + xp) + 2) = 0.0f;
					}
				}
			}
		}
		address++;
	}
}


void Quilting::matchFeatures (int currentIteration) {
	flann::KDTreeIndexParams kdIndexParams (params.trees);
	flann::Index kdIndex (features, kdIndexParams);
	int step = (currentIteration == 0 ? 2 * radius - overlap + 1 : radius + 1);
	int shift = 3 * (2 * radius + 1) * (2 * radius + 1);
	int shift_2 = 3 * overlap * (2 * radius + 1);
	for (int i = radius ; i <= targetNormalMap.rows - 1 - radius ; i += step) {
		for (int j = radius ; j <= targetNormalMap.cols - 1 - radius ; j += step) {
			vector<int> kdCurrentIndex (1, -1);
		    vector<float> kdCurrentDist (1, FLT_MAX);
		    vector<float> query (dimensionality, 0.0f);
			/* Filling in the query vector, first with the normals */
			for (int y = max(0, i - radius) ; y <= min(targetNormalMap.rows - 1, i + radius) ; y++) {
				for (int x = max(0, j - radius) ; x <= min(targetNormalMap.cols - 1, j + radius) ; x++) {
					int yp = y - (i - radius);
					int xp = x - (j - radius);
					query[3 * (yp * ((radius << 1) + 1) + xp) + 0] = sqrt(alpha) * targetNormalMap.at<Vec4f>(y, x)[0];
					query[3 * (yp * ((radius << 1) + 1) + xp) + 1] = sqrt(alpha) * targetNormalMap.at<Vec4f>(y, x)[1];
					query[3 * (yp * ((radius << 1) + 1) + xp) + 2] = sqrt(alpha) * targetNormalMap.at<Vec4f>(y, x)[2];
				}
			}
			/* Then with the colors */
			if (currentIteration == 0) {
				for (int y = max(0, i - radius) ; y <= min(targetNormalMap.rows - 1, i + radius) ; y++) {
					for (int x = max(0, j - radius) ; x <= min(targetNormalMap.cols - 1, j + radius) ; x++) {
						int yp = y - (i - radius);
						int xp = x - (j - radius);
						if (yp < overlap) {
							query[shift + 3 * (yp * ((radius << 1) + 1) + xp) + 0] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[0];
							query[shift + 3 * (yp * ((radius << 1) + 1) + xp) + 1] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[1];
							query[shift + 3 * (yp * ((radius << 1) + 1) + xp) + 2] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[2];
						} else if (xp < overlap) {
							query[shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 0] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[0];
							query[shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 1] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[1];
							query[shift + shift_2 + 3 * ((yp - overlap) * overlap + xp) + 2] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[2];
						}
					}
				}
			} else {
				for (int y = max(0, i - radius) ; y <= min(targetNormalMap.rows - 1, i + radius) ; y++) {
					for (int x = max(0, j - radius) ; x <= min(targetNormalMap.cols - 1, j + radius) ; x++) {
						int yp = y - (i - radius);
						int xp = x - (j - radius);
						query[shift + 3 * (yp * ((radius << 1) + 1) + xp) + 0] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[0];
						query[shift + 3 * (yp * ((radius << 1) + 1) + xp) + 1] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[1];
						query[shift + 3 * (yp * ((radius << 1) + 1) + xp) + 2] = sqrt(beta) * newAlbedo.at<Vec4f>(y, x)[2];
					}
				}
			}
			/* K-d tree search */
			kdIndex.knnSearch (query, kdCurrentIndex, kdCurrentDist, 1, flann::SearchParams(params.recursions));
			/* Paste */
			// ...
		}
	}	
}

