#ifndef _QUILTING_H_
#define _QUILTING_H_

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <deque>
#include <utility>
#include "method.h"

using namespace cv;
using namespace std;

enum PasteMode {
	NONE,
	TOP,
	LEFT,
	CORNER,
	ENTIRE
};

class Quilting : public Method
{
public:
	Quilting (Mat & _sourceNormalMap,
			  Mat & _sourceAlbedo,
			  Mat & _sourceShading,
			  Mat & _targetNormalMap,
			  Mat & _targetShading,
			  Mat & _newAlbedo,
			  Mat & _newShading,
			  Mat & _newImage,
			  MethodParams & _params);

	~Quilting ();

	void textureTransfer ();

protected:
	int detectFeatures ();

	void buildFeatures (int currentIteration);

	void matchFeatures (int currentIteration);

	void paste (int index, Point2i & center, enum PasteMode mode);

	void verticalBoundaryCut (int index, Point2i & target, vector<int> & coords);

	void horizontalBoundaryCut (int index, Point2i & target, vector<int> & coords);

	int radius;
	int overlap;
	int entries;
	int dimensionality;
	float alpha;
	float beta;
	deque<Point2i> centers;
	Mat features;

	Mat diagram;
};

#endif
