#ifndef _METHOD_H_
#define _METHOD_H_

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "arg.h"

using namespace cv;

class Method {
public:
	Method (Mat & _sourceNormalMap,
			Mat & _sourceAlbedo,
			Mat & _sourceShading,
			Mat & _targetNormalMap,
			Mat & _targetShading,
			Mat & _newAlbedo,
			Mat & _newShading,
			MethodParams & _params);

	~Method ();
	
	virtual void textureTransfer () = 0;

protected:
	Mat sourceNormalMap;
	Mat sourceAlbedo;
	Mat sourceShading;
	Mat targetNormalMap;
	Mat targetShading;
	Mat newAlbedo;
	Mat newShading;
	MethodParams params;
};

#endif
