#ifndef _QUILTING_H_
#define _QUILTING_H_

#include "method.h"

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
	void buildFeatures ();

	void matchFeatures ();
};

#endif
