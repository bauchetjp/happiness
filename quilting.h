#ifndef _QUILTING_H_
#define _QUILTING_H_

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <list>
#include <utility>
#include "method.h"

using namespace cv;
using namespace std;

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
			  MethodParams & _params);

	~Quilting ();

	void textureTransfer ();

protected:
	int detectFeatures ();

	void buildFeatures (int currentIteration);

	void matchFeatures (int currentIteration);

	int radius;
	int overlap;
	int entries;
	int dimensionality;
	float alpha;
	float beta;
	list<Point2i> centers;
	Mat features;
};

#endif
