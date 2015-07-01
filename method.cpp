#include "method.h"

Method::Method (Mat & _sourceNormalMap,
				Mat & _sourceAlbedo,
				Mat & _sourceShading,
				Mat & _targetNormalMap,
				Mat & _targetShading,
				Mat & _newAlbedo,
				Mat & _newShading,
				MethodParams & _params) :
	sourceNormalMap (_sourceNormalMap),
	sourceAlbedo (_sourceAlbedo),
	sourceShading (_sourceShading),
	targetNormalMap (_targetNormalMap),
	targetShading (_targetShading),
	newAlbedo (_newAlbedo),
	newShading (_newShading),
	params (_params)
{

}

Method::~Method () {

}

