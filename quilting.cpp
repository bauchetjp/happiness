#include "method.h"
#include "quilting.h"

using namespace cv;

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

}


Quilting::~Quilting() {

}


void Quilting::textureTransfer () {

}


void Quilting::buildFeatures () {

}


void Quilting::matchFeatures () {

}
