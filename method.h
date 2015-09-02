#ifndef _METHOD_H_
#define _METHOD_H_

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <vector>

#include "arg.h"

using namespace cv;
using namespace std;

class Method {
public:
	Method (Mat & _sourceNormalMap,
			Mat & _sourceAlbedo,
			Mat & _sourceShading,
			Mat & _targetNormalMap,
			Mat & _targetShading,
			Mat & _newAlbedo,
			Mat & _newShading,
			Mat & _newImage,
			MethodParams & _params);

	void save ();

	virtual ~Method ();
	
	virtual void textureTransfer () = 0;

	void shadingTransfer ();

protected:
	Mat sourceNormalMap;
	Mat sourceAlbedo;
	Mat sourceShading;
	Mat targetNormalMap;
	Mat targetShading;
	Mat newAlbedo;
	Mat newShading;
	Mat newImage;
	MethodParams params;

	void _save (string filename, Mat m);

private:
	void updateParameters(Mat & frame, vector<double> & stats);

	void computeSaturationAngle(Mat & shadingMap, Mat & normalMap, float & saturationAngle, bool verbose, float threshold = 0.8);

	void linearHistogramRemapping(float slope);

	void nonLinearHistogramRemapping(float slope, float abscissa, float exponent);

	void generateOptimalShadingsPos(float slope, float absStep, float expMin, float expStep, float expMax, string prefix);

	void generateOptimalShadingsNeg(float slope, float absStep, float expMin, float expStep, float expMax, string prefix);

	void multiply(Mat & albedo, Mat & shading, Mat & result);

	void context(FILE* logFile);

	float tab[6];

	vector<double> stats_shading_s;
	vector<double> stats_shading_t;
	vector<double> stats_shading_n;
	float angle_saturation_s;
	float angle_saturation_n;
};

float dot(Vec3f & u, Vec3f & v); 

float length(Vec3f & u); 

#endif
