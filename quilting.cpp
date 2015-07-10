#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/flann/flann.hpp>
#include <deque>
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
	centers = deque<Point2i>();
}


Quilting::~Quilting() {
	centers.clear();
	features.release();
}


void Quilting::textureTransfer () {
	for (int i = 0 ; i < params.quiltingIterations ; i++) {
		if (i == 0) {
			radius = params.quiltingRadius;
		} else {
			radius = max(2, (2 * radius) / 3);
		}
		overlap = max(2, radius / 3);	
		if (i == 0) {		
			dimensionality = 3 * ((2 * radius + 1) * (2 * radius + 1) + 2 * overlap * (2 * radius + 1) - overlap * overlap);
		} else {
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
	deque<Point2i>::iterator it;
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
	int requests = 0;
	for (int i = radius ; i <= targetNormalMap.rows - 1 - radius ; i += step) {
		cout << currentIteration << " " << i << endl;
		for (int j = radius ; j <= targetNormalMap.cols - 1 - radius ; j += step) {
			requests++;
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
			Point2i center (j, i);
			enum PasteMode mode;
			if (currentIteration == 0) {
				if (i == radius && j == radius) {
					mode = NONE;
				} else if (i == radius) {
					mode = LEFT;
				} else if (j == radius) {
					mode = TOP;
				} else {
					mode = CORNER;
				}
			} else {
				mode = ENTIRE;
			}
			paste (kdCurrentIndex[0], center, mode);			
		}
	}	
	cout << requests << endl;
}


void Quilting::paste (int index, Point2i & target, enum PasteMode mode) {
	Point2i source = centers[index];
	if (mode == NONE) {
		for (int i = -radius ; i <= radius ; i++) {
			for (int j = -radius ; j <= radius ; j++) {
				int t_i = clamp(0, target.y + i, targetNormalMap.rows - 1);
				int t_j = clamp(0, target.x + j, targetNormalMap.cols - 1);
				int s_i = clamp(0, source.y + i, sourceNormalMap.rows - 1);
				int s_j = clamp(0, source.x + j, sourceNormalMap.cols - 1);
				if (targetNormalMap.at<Vec4f>(t_i, t_j)[3] > 0.9f) {
					newAlbedo.at<Vec4f>(t_i, t_j) = sourceAlbedo.at<Vec4f>(s_i, s_j);
					newAlbedo.at<Vec4f>(t_i, t_j)[3] = 1.0f;	
				} else {
					newAlbedo.at<Vec4f>(t_i, t_j) = Vec4f(0.0f, 0.0f, 0.0f, 0.0f);
				}
			}
		}
	} else if (mode == TOP) {
		vector<int> coords (2 * radius + 1, -1);
		horizontalBoundaryCut (index, target, coords);
		for (int j = -radius ; j <= radius ; j++) {
			for (int i = -radius + coords[radius + j] ; i <= radius ; i++) {
				int t_i = clamp(0, target.y + i, targetNormalMap.rows - 1);
				int t_j = clamp(0, target.x + j, targetNormalMap.cols - 1);
				int s_i = clamp(0, source.y + i, sourceNormalMap.rows - 1);
				int s_j = clamp(0, source.x + j, sourceNormalMap.cols - 1);
				if (targetNormalMap.at<Vec4f>(t_i, t_j)[3] > 0.9f) {
					newAlbedo.at<Vec4f>(t_i, t_j) = sourceAlbedo.at<Vec4f>(s_i, s_j);
					newAlbedo.at<Vec4f>(t_i, t_j)[3] = 1.0f;	
				} else {
					newAlbedo.at<Vec4f>(t_i, t_j) = Vec4f(0.0f, 0.0f, 0.0f, 0.0f);
				}
			}
		}
	} else if (mode == LEFT) {
		vector<int> coords (2 * radius + 1, -1);
		verticalBoundaryCut (index, target, coords);
		for (int i = -radius ; i <= radius ; i++) {
			for (int j = -radius + coords[radius + i] ; j <= radius ; j++) {
				int t_i = clamp(0, target.y + i, targetNormalMap.rows - 1);
				int t_j = clamp(0, target.x + j, targetNormalMap.cols - 1);
				int s_i = clamp(0, source.y + i, sourceNormalMap.rows - 1);
				int s_j = clamp(0, source.x + j, sourceNormalMap.cols - 1);
				if (targetNormalMap.at<Vec4f>(t_i, t_j)[3] > 0.9f) {
					newAlbedo.at<Vec4f>(t_i, t_j) = sourceAlbedo.at<Vec4f>(s_i, s_j);
					newAlbedo.at<Vec4f>(t_i, t_j)[3] = 1.0f;	
				} else {
					newAlbedo.at<Vec4f>(t_i, t_j) = Vec4f(0.0f, 0.0f, 0.0f, 0.0f);
				}
			}
		}
	} else if (mode == CORNER || mode == ENTIRE) {
		vector<int> hcoords (2 * radius + 1, -1);
		vector<int> vcoords (2 * radius + 1, -1);
		horizontalBoundaryCut (index, target, hcoords);
		verticalBoundaryCut (index, target, vcoords);
		for (int i = -radius ; i <= radius ; i++) {
			for (int j = -radius ; j <= radius ; j++) {
				if (i + radius >= hcoords[j + radius] && j + radius >= vcoords[i + radius]) {
					int t_i = clamp(0, target.y + i, targetNormalMap.rows - 1);
					int t_j = clamp(0, target.x + j, targetNormalMap.cols - 1);
					int s_i = clamp(0, source.y + i, sourceNormalMap.rows - 1);
					int s_j = clamp(0, source.x + j, sourceNormalMap.cols - 1);
					if (targetNormalMap.at<Vec4f>(t_i, t_j)[3] > 0.9f) {
						newAlbedo.at<Vec4f>(t_i, t_j) = sourceAlbedo.at<Vec4f>(s_i, s_j);
						newAlbedo.at<Vec4f>(t_i, t_j)[3] = 1.0f;	
					} else {
						newAlbedo.at<Vec4f>(t_i, t_j) = Vec4f(0.0f, 0.0f, 0.0f, 0.0f);
					}
				}
			}
		}
	}
}


void Quilting::horizontalBoundaryCut (int index, Point2i & target, vector<int> & coords) {
	Mat e (overlap, 2 * radius + 1, CV_32F);
	Mat E (overlap, 2 * radius + 1, CV_32F);
	Mat selected (1, dimensionality, CV_32F);
	features.row(index).copyTo(selected.row(0));
	int shift = 3 * (2 * radius + 1) * (2 * radius + 1);
	for (int q = 0 ; q < e.rows ; q++) {
		for (int p = 0 ; p < e.cols ; p++) {
			Vec3f u, v;
			int y = clamp(0, target.y - radius + q, newAlbedo.rows - 1);
			int x = clamp(0, target.x - radius + p, newAlbedo.cols - 1);
			for (int d = 0 ; d < 3 ; d++) {
				u[d] = newAlbedo.at<Vec4f>(y, x)[d];
				v[d] = features.at<float>(0, shift + 3 * (q * (2 * radius + 1) + p) + d);
			}
			e.at<float>(q, p) = 0.0;
			for (int d = 0 ; d < 3 ; d++) {
				e.at<float>(q, p) = sqrt((u[d] - v[d]) * (u[d] - v[d]));
			}
		}
	} 
	for (int q = 0 ; q < e.rows ; q++) {
		E.at<float>(q, 0) = e.at<float>(q, 0);
	}
	for (int p = 1 ; p < e.cols ; p++) {
		for (int q = 0 ; q < e.rows ; q++) {
			if (q == 0) {
				E.at<float>(0, p) = e.at<float>(0, p) + min(E.at<float>(0, p - 1), E.at<float>(1, p - 1));
			} else if (q == overlap - 1) {
				E.at<float>(overlap - 1, p) = e.at<float>(overlap - 1, p) + min(E.at<float>(overlap - 2, p - 1), E.at<float>(overlap - 1, p - 1));
			} else {
				E.at<float>(q, p) = e.at<float>(q, p) + min(E.at<float>(q, p - 1), min(E.at<float>(q - 1, p - 1), E.at<float>(q + 1, p - 1)));
			}
		}
	}
	float columnMin = FLT_MAX;
	for (int q = 0 ; q < e.rows ; q++) {
		if (E.at<float>(q, e.cols - 1) < columnMin) {
			coords[e.cols - 1] = q;
			columnMin = E.at<float>(q, e.cols - 1);
		}
	}
	for (int p = e.cols - 2 ; p >= 0 ; p--) {
		float expectedMin = E.at<float>(coords[p + 1], p + 1) - e.at<float>(coords[p + 1], p + 1);
		for (int q = max(0, coords[p + 1] - 1) ; q <= min(overlap - 1, coords[p + 1] + 1) ; q++) {
			if (fabs(E.at<float>(q, p) - expectedMin) < 1e-5) {
				coords[p] = q;
				break;
			}
		}
	}
}


void Quilting::verticalBoundaryCut (int index, Point2i & target, vector<int> & coords) {
	Mat e (2 * radius + 1, overlap, CV_32F);
	Mat E (2 * radius + 1, overlap, CV_32F);
	Mat selected (1, dimensionality, CV_32F);
	features.row(index).copyTo(selected.row(0));
	int shift = 3 * (2 * radius + 1) * (2 * radius + 1);
	int shift_2 = 3 * overlap * (2 * radius + 1);
	for (int q = 0 ; q < e.rows ; q++) {
		for (int p = 0 ; p < e.cols ; p++) {
			Vec3f u, v;
			int y = clamp(0, target.y - radius + q, newAlbedo.rows - 1);
			int x = clamp(0, target.x - radius + p, newAlbedo.cols - 1);
			for (int d = 0 ; d < 3 ; d++) {
				u[d] = newAlbedo.at<Vec4f>(y, x)[d];
				if (q < overlap) {
					v[d] = features.at<float>(0, shift + 3 * (q * ((radius << 1) + 1) + p) + d);
				} else {
					v[d] = features.at<float>(0, shift + shift_2 + 3 * (((q - overlap) * overlap) + p) + d);
				}
			}
			e.at<float>(q, p) = 0.0;
			for (int d = 0 ; d < 3 ; d++) {
				e.at<float>(q, p) = sqrt((u[d] - v[d]) * (u[d] - v[d]));
			}
		}
	}
	for (int p = 0 ; p < e.cols ; p++) {
		E.at<float>(0, p) = e.at<float>(0, p);
	}
	for (int q = 1 ; q < e.rows ; q++) {
		for (int p = 0 ; p < e.cols ; p++) {
			if (p == 0) {
				E.at<float>(q, 0) = e.at<float>(q, 0) + min(E.at<float>(q - 1, 0), E.at<float>(q - 1, 1));
			} else if (p == e.cols - 1) {
				E.at<float>(q, e.cols - 1) = e.at<float>(q, e.cols - 1) + min(E.at<float>(q - 1, e.cols - 2), E.at<float>(q - 1, e.cols - 1));
			} else {
				E.at<float>(q, p) = e.at<float>(q, p) + min(E.at<float>(q - 1, p), min(E.at<float>(q - 1, p - 1), E.at<float>(q - 1, p + 1)));
			}
		}
	}
	float rowMin = FLT_MAX;
	for (int p = 0 ; p < e.cols ; p++) {
		if (E.at<float>(e.rows - 1, p) < rowMin) {
			coords[e.rows - 1] = p;
			rowMin = E.at<float>(e.rows - 1, p);
		}
	}
	for (int q = e.rows - 2 ; q >= 0 ; q--) {
		float expectedMin = fabs(E.at<float>(q + 1, coords[q + 1]) - e.at<float>(q + 1, coords[q + 1]));
		for (int p = max(0, coords[q + 1] - 1) ; p < min(e.rows - 1, coords[q + 1] + 1) ; p++) {
			if (fabs(E.at<float>(q, p) - expectedMin) < 1e-5) {
				coords[q] = p;
				break;
			}
		}
	}
}

