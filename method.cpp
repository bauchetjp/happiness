#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "method.h"

#define minimum(a, b) ((a) <= (b) ? (a) : (b))
#define maximum(a, b) ((a) >= (b) ? (a) : (b))
#define clamp(a, x, b) ((x) <= (a) ? (a) : ((b) <= (x) ? (b) : (x)))


Method::Method (Mat & _sourceNormalMap,
				Mat & _sourceAlbedo,
				Mat & _sourceShading,
				Mat & _targetNormalMap,
				Mat & _targetShading,
				Mat & _newAlbedo,
				Mat & _newShading,
				Mat & _newImage,
				MethodParams & _params) :
	sourceNormalMap (_sourceNormalMap),
	sourceAlbedo (_sourceAlbedo),
	sourceShading (_sourceShading),
	targetNormalMap (_targetNormalMap),
	targetShading (_targetShading),
	newAlbedo (_newAlbedo),
	newShading (_newShading),
	newImage (_newImage),
	params (_params)
{

}


Method::~Method () {

}


void Method::save () {
	_save(params.prefix + "_sourceNormalMap.png", sourceNormalMap);
	_save(params.prefix + "_sourceAlbedo.png", sourceAlbedo);
	_save(params.prefix + "_targetNormalMap.png", targetNormalMap);
	_save(params.prefix + "_target.png", newAlbedo);
}


void Method::_save (string filename, Mat m) {
	Mat x (480, 640, CV_8UC4);
	for (int i = 0 ; i < 480 ; i++) {
		for (int j = 0 ; j < 640 ; j++) {
			Vec4b & bgra = x.at<Vec4b>(i, j);
			bgra[0] = int(clamp(0, 255 * m.at<Vec4f>(i, j)[2], 255));
			bgra[1] = int(clamp(0, 255 * m.at<Vec4f>(i, j)[1], 255));
			bgra[2] = int(clamp(0, 255 * m.at<Vec4f>(i, j)[0], 255));
			bgra[3] = int(clamp(0, 255 * m.at<Vec4f>(i, j)[3], 255));
		}
	}
	imwrite(filename, x);
}


float dot(Vec3f & u, Vec3f & v) {
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}


float length(Vec3f & u) {
	return sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
}


void Method::shadingTransfer() {
	updateParameters (sourceShading, stats_shading_s);
	updateParameters (targetShading, stats_shading_t);
	computeSaturationAngle (sourceShading, sourceNormalMap, angle_saturation_s, false);
	if (!params.automaticShading) {
		nonLinearHistogramRemapping (params.alpha, params.beta, params.gamma);
		multiply(newAlbedo, newShading, newImage);
		_save(params.prefix + "_newShading.png", newShading);
		_save(params.prefix + "_result.png", newImage);
	} else {
		linearHistogramRemapping (params.alpha);
		if (stats_shading_s[3] > stats_shading_n[3]) {
			generateOptimalShadingsPos (params.alpha, 0.025f, 0.1f, 0.1f, 100.0f, params.prefix);
		} else {
			generateOptimalShadingsNeg (params.alpha, 0.025f, -100.0f, 0.1f, -0.1f, params.prefix);
		}
	}
}


void Method::updateParameters (Mat & frame, vector<float> & stats) {
	stats.clear();
	int n = 0;
    float min = FLT_MAX;
    float max = -FLT_MAX;
    float mean = 0;
    float sd = 0;
    float skew = 0;
    for (int y = 0 ; y < frame.rows ; y++) {
        for (int x = 0 ; x < frame.cols ; x++) {
            if (frame.at<Vec4f>(y, x)[3] > 1.0e-6) {
				float r = frame.at<Vec4f>(y, x)[0];
                mean += r;
                if (r > max) {
                    max = r;
                } else if (r < min) {
                    min = r;
                }
                n++;
            }
        }
    }
    mean /= n;
    for (int y = 0 ; y < frame.rows ; y++) {
        for (int x = 0 ; x < frame.cols ; x++) {
            if (frame.at<Vec4f>(y, x)[3] > 1.0e-6) {
                float r = frame.at<Vec4f>(y, x)[0];
                sd += (r - mean) * (r - mean);
            }
        }
    }
    sd /= n;
    sd = sqrt(sd);
    for (int y = 0 ; y < frame.rows ; y++) {
        for (int x = 0 ; x < frame.cols ; x++) {
            if (frame.at<Vec4f>(y, x)[3] > 1.0e-6) {
                float r = frame.at<Vec4f>(y, x)[0];
                skew += pow((r - mean) / sd, 3.0);
            }
        }
    }
    skew /= float(n);
    stats.push_back(float(n));    // 0 : Realizations
    stats.push_back(mean);        // 1 : Mean
    stats.push_back(sd);          // 2 : Standard deviation
    stats.push_back(skew);        // 3 : Skewness
    stats.push_back(min);         // 4 : Minimum
	stats.push_back(max);         // 5 : Maximum
}


void Method::computeSaturationAngle(Mat & shadingMap, Mat & normalMap, float & saturationAngle, bool verbose, float threshold) {
	// Finding the highest value in the layer and its argmax
    float max_quantity = -FLT_MAX;
    std::list<Vec4f> argmaxList;
	int height = shadingMap.rows;
	int width = shadingMap.cols;
    for (int y = 0 ; y < height ; y++) {
        for (int x = 0 ; x < width ; x++) {
            if (shadingMap.at<Vec4f>(y, x)[3] > 1.0e-6) {
                float quantity = shadingMap.at<Vec4f>(y, x)[0];
                if (quantity > max_quantity) {
                    max_quantity = quantity;
                    argmaxList.clear();
                    argmaxList.push_back(normalMap.at<Vec4f>(y, x));
                } else if (fabs(quantity - max_quantity) < 1.0e-4) {
                    argmaxList.push_back(normalMap.at<Vec4f>(y, x));
                }
            }
        }
    }
	// In the case there are several maximums, compute the average normal
    Vec3f argmax (0.0, 0.0, 0.0);
    std::list<Vec4f>::iterator it1;
    for (it1 = argmaxList.begin() ; it1 != argmaxList.end() ; it1++) {
        argmax += Vec3f((*it1)[0], (*it1)[1], (*it1)[2]);
    }
    argmax = (1.0f / float(argmaxList.size())) * argmax;
	// Retrieve the areas of the image where the value in the layer is greater than a certain threshold
    // and return the borders of these areas, using morphologic operators
    // Initialization
    float boundary = minimum(1.0, threshold * max_quantity);
	Mat borders (height, width, CV_32S);
	Mat erodedBorders (height, width, CV_32S);
	Mat erodedBordersInt (height, width, CV_32S);
	for (int y = 0 ; y < height ; y++) {
		for (int x = 0 ; x < width ; x++) {
			if (shadingMap.at<Vec4f>(y, x)[0] < boundary) {
                borders.at<int>(y, x) = erodedBorders.at<int>(y, x) = erodedBordersInt.at<int>(y, x) = 0;
            } else {
                borders.at<int>(y, x) = erodedBorders.at<int>(y, x) = erodedBordersInt.at<int>(y, x) = 1;
            }
		}
	}
	// Erosion
	int size = 1;
    for (int i = size ; i < height - size ; i++) {
        for (int j = size ; j < width - size ; j++) {
            int m = borders.at<int>(i, j);
            for (int k = i - size ; k <= i + size ; k++) {
                if (borders.at<int>(k, j) < m) {
                    m = borders.at<int>(k, j);
                }
            }
            erodedBordersInt.at<int>(i, j) = m;
        }
    }
    for (int j = size ; j < width - size ; j++) {
        for (int i = size ; i < height - size ; i++) {
            int m = erodedBordersInt.at<int>(i, j);
            for (int l = j - size ; l <= j + size ; l++) {
                if (erodedBordersInt.at<int>(i, l) < m) {
                    m = erodedBordersInt.at<int>(i, l);
                }
            }
            erodedBorders.at<int>(i, j) = m;
        }
    }
    //	erode (borders, erodedBorders, Mat());
	// Compute the borders of the areas of high shading
    for (int y = 0 ; y < height ; y++) {
        for (int x = 0 ; x < width ; x++) {
            borders.at<int>(y, x) -= erodedBorders.at<int>(y, x); // borders = 0 or 1
        }
    }
	// Now, loop on every entry of the normal map located on the border
    // and calculate the average angle between the read normal and the reference normal
    float angle = 0.0;
    int terms = 0;
    for (int y = 0 ; y < height ; y++) {
        for (int x = 0 ; x < width ; x++) {
            if (borders.at<int>(y, x)) {
                Vec3f normal = Vec3f(normalMap.at<Vec4f>(y, x)[0], normalMap.at<Vec4f>(y, x)[1], normalMap.at<Vec4f>(y, x)[2]);
                angle += acos(dot(argmax, normal) / (length(argmax) * length(normal)));
                terms++;
            }
        }
    }
    angle /= terms;
    saturationAngle = angle;
}


void Method::linearHistogramRemapping(float slope) {
    // Applies the transformation
    // The target histogram is brought into correspondance
    for (int y = 0 ; y <= targetShading.rows - 1 ; y++) {
        for (int x = 0 ; x <= targetShading.cols - 1 ; x++) {
            if (targetShading.at<Vec4f>(y, x)[3] > 1.0e-6) {
                // mean_shading_s + slope * (sd_shading_s / sd_shading_t) * (x - mean_shading_t);
                newShading.at<Vec4f>(y, x)[0] = stats_shading_s[1] + slope * (stats_shading_s[2] / stats_shading_t[2]) * (targetShading.at<Vec4f>(y, x)[0] - stats_shading_t[1]);
                newShading.at<Vec4f>(y, x)[1] = stats_shading_s[1] + slope * (stats_shading_s[2] / stats_shading_t[2]) * (targetShading.at<Vec4f>(y, x)[1] - stats_shading_t[1]);
                newShading.at<Vec4f>(y, x)[2] = stats_shading_s[1] + slope * (stats_shading_s[2] / stats_shading_t[2]) * (targetShading.at<Vec4f>(y, x)[2] - stats_shading_t[1]);
                newShading.at<Vec4f>(y, x)[3] = 1.0;
				
            } else {
                newShading.at<Vec4f>(y, x) = Vec4f(0.0, 0.0, 0.0, 0.0);
            }
        }
    }
    updateParameters(newShading, stats_shading_n);
    computeSaturationAngle(newShading, targetNormalMap, angle_saturation_n, false);
}


void Method::nonLinearHistogramRemapping(float slope, float abscissa, float exponent) {
    float break_t;
    if (abscissa >= 0) {
        break_t = (1 - abscissa) * stats_shading_t[1] + abscissa * stats_shading_t[5]; // max
    } else {
        break_t = (1 + abscissa) * stats_shading_t[1] - abscissa * stats_shading_t[4]; // min
    }
    for (int y = 0 ; y <= targetShading.rows - 1 ; y++) {
        for (int x = 0 ; x <= targetShading.cols - 1 ; x++) {
            if (targetShading.at<Vec4f>(y, x)[3] > 1.0e-6) {
                if (targetShading.at<Vec4f>(y, x)[0] <= break_t) {
                    // Linear portion
                    newShading.at<Vec4f>(y, x)[0] = stats_shading_s[1] + slope * (stats_shading_s[2] / stats_shading_t[2]) * (targetShading.at<Vec4f>(y, x)[0] - stats_shading_t[1]);
                    newShading.at<Vec4f>(y, x)[1] = stats_shading_s[1] + slope * (stats_shading_s[2] / stats_shading_t[2]) * (targetShading.at<Vec4f>(y, x)[1] - stats_shading_t[1]);
                    newShading.at<Vec4f>(y, x)[2] = stats_shading_s[1] + slope * (stats_shading_s[2] / stats_shading_t[2]) * (targetShading.at<Vec4f>(y, x)[2] - stats_shading_t[1]);
                    newShading.at<Vec4f>(y, x)[3] = 1.0;
                } else {
                    // Exponential portion
                    float y_coord = stats_shading_s[1] + slope * (stats_shading_s[2] / stats_shading_t[2]) * (break_t - stats_shading_t[1] - 1.0 / exponent);
                    float f = ((slope * stats_shading_s[2]) / (exponent * stats_shading_t[2]));
                    newShading.at<Vec4f>(y, x)[0] = y_coord + f * exp(exponent * (targetShading.at<Vec4f>(y, x)[0] - break_t));
                    newShading.at<Vec4f>(y, x)[1] = y_coord + f * exp(exponent * (targetShading.at<Vec4f>(y, x)[1] - break_t));
                    newShading.at<Vec4f>(y, x)[2] = y_coord + f * exp(exponent * (targetShading.at<Vec4f>(y, x)[2] - break_t));
                    newShading.at<Vec4f>(y, x)[3] = 1.0;
                }
            } else {
                newShading.at<Vec4f>(y, x) = Vec4f(0.0, 0.0, 0.0, 0.0);
            }
        }
    }
    updateParameters(newShading, stats_shading_n);
    computeSaturationAngle(newShading, targetNormalMap, angle_saturation_n, false);
}


void Method::context(FILE* logFile) {
    fprintf(logFile, "Source shading : \n");
    fprintf(logFile, "    n_s   = %i\n", int(stats_shading_s[0]));
    fprintf(logFile, "    min   = %f\n", stats_shading_s[4]);
    fprintf(logFile, "    max   = %f\n", stats_shading_s[5]);
    fprintf(logFile, "    mean  = %f\n", stats_shading_s[1]);
    fprintf(logFile, "    sd    = %f\n", stats_shading_s[2]);
    fprintf(logFile, "    skew  = %f\n", stats_shading_s[3]);
    fprintf(logFile, "    angle = %f\n\n", angle_saturation_s);
    fprintf(logFile, "Target shading : \n");
    fprintf(logFile, "    n_t   = %i\n", int(stats_shading_s[0]));
    fprintf(logFile, "    min   = %f\n", stats_shading_s[4]);
    fprintf(logFile, "    max   = %f\n", stats_shading_s[5]);
    fprintf(logFile, "    mean  = %f\n", stats_shading_t[1]);
    fprintf(logFile, "    sd    = %f\n", stats_shading_t[2]);
    fprintf(logFile, "    skew  = %f\n", stats_shading_t[3]);
}


void Method::generateOptimalShadingsPos(float slope, float absStep, float expMin, float expStep, float expMax, string prefix) {
	int expItMax = int(fabs(expMax - expMin) / expStep) - 1;
    int expArg = 0;
    int absIterations = int(2.0 / absStep);
   	string logFilename (prefix);
	logFilename.append(".log");
    float minDiffAnglesSh = FLT_MAX;
    int argminDiffAnglesSh = -1;
    FILE* logFile = fopen (logFilename.c_str(), "w");
    if (logFile != NULL) {
        fprintf(logFile, "** Optimizes the skewness of the shading layer... **\n\n");
        context(logFile);
        for (int absIt = 0 ; absIt < absIterations ; absIt++) {
            // We would like to find an isocurve
            float abs = -1.0 + absIt * absStep;
            int expIt = expArg;
            float expCur  = expMin + expIt * expStep;
            float expNext = expMin + (expIt + 1) * expStep;
            nonLinearHistogramRemapping(slope, abs, expNext);
            while (stats_shading_n[3] <= stats_shading_s[3] && expIt < expItMax - 1) {
                expIt++;
                expCur = expNext;
                expNext = expMin + (expIt + 1) * expStep;
                nonLinearHistogramRemapping(slope, abs, expNext);
            }
            nonLinearHistogramRemapping(slope, abs, expCur);
            float distCur = fabs(stats_shading_n[3] - stats_shading_s[3]);
            bool sup = (stats_shading_n[3] <= stats_shading_s[3]);
            nonLinearHistogramRemapping(slope, abs, expNext);
            float distNext = fabs(stats_shading_n[3] - stats_shading_s[3]);
            bool inf = (stats_shading_s[3] < stats_shading_n[3]);
			float expRes;
            if (distCur < distNext) {
                expRes = expCur;
            } else {
                expRes = expNext;
            }
            expArg = expIt;
            float diffAnglesSh = fabs(angle_saturation_s - angle_saturation_n);
            if (diffAnglesSh < minDiffAnglesSh) {
                minDiffAnglesSh = diffAnglesSh;
                argminDiffAnglesSh = absIt;
            }
            multiply(newAlbedo, newShading, newImage);
			string imageFilename (prefix), shadingFilename(prefix);
			imageFilename += "_res_"; shadingFilename += "_shading_";
			imageFilename += to_string(absIt); shadingFilename += to_string(absIt);
			imageFilename += "_"; shadingFilename += "_";
			imageFilename += to_string(int(100 * slope)); shadingFilename += to_string(int(100 * slope));
			imageFilename += "_"; shadingFilename += "_";
			imageFilename += to_string(int(100 * abs)); shadingFilename += to_string(int(100 * abs));
			imageFilename += "_"; shadingFilename += "_";
			imageFilename += to_string(int(100 * expRes)); shadingFilename += to_string(int(100 * expRes));
			imageFilename += ".png"; shadingFilename += ".png";
            _save(imageFilename, newImage);
			_save(shadingFilename, newShading);
            fprintf(logFile, "It. %2i (break = %6.3f) : Minimum found in exp = %5.1f : ", absIt, abs, expRes);
            fprintf(logFile, "mean = %9.5f, sd = %9.5f, skew = %9.5f, angle = %9.5f",
                    stats_shading_n[1], stats_shading_n[2], stats_shading_n[3], angle_saturation_n);
            if (!inf || !sup) {
                fprintf(logFile, " *");
            }
            fprintf(logFile, "\n");
        }
        fprintf(logFile, "\n");
        fprintf(logFile, "Optimal suggested result using shading saturation : %i\n", argminDiffAnglesSh);
        fclose(logFile);
    } else {
        std::cout << "generateOptimalShadingsPos() : error : couldn't open a log file" << std::endl;
    }
}


void Method::generateOptimalShadingsNeg(float slope, float absStep, float expMin, float expStep, float expMax, string prefix) {
int expItMax = int(fabs(expMax - expMin) / expStep) - 1;
    int expArg = expItMax - 1;
    int absIterations = int(2.0 / absStep);
    string logFilename (prefix);
	logFilename += ".log";
    float minDiffAnglesSh = FLT_MAX;
    int argminDiffAnglesSh = -1;
    FILE* logFile = fopen (logFilename.c_str(), "w");
    if (logFile != NULL) {
        fprintf(logFile, "** Optimizes the skewness of the shading layer... **\n\n");
        context(logFile);
        for (int absIt = 0 ; absIt <= absIterations ; absIt++) {
            float abs = -1.0 + absIt * absStep;
            int expIt = expArg;
            float expCur  = expMin + expIt * expStep;
            float expNext = expMin + (expIt - 1) * expStep;
            nonLinearHistogramRemapping(slope, abs, expNext);
            while (stats_shading_n[3] >= stats_shading_s[3] && expIt > 0) {
                expIt--;
                expCur = expNext;
                expNext = expMin + (expIt - 1) * expStep;
                nonLinearHistogramRemapping(slope, abs, expNext);
            }
            nonLinearHistogramRemapping(slope, abs, expCur);
            float distCur = fabs(stats_shading_n[3] - stats_shading_s[3]);
            bool sup = (stats_shading_n[3] >= stats_shading_s[3]);
            nonLinearHistogramRemapping(slope, abs, expNext);
            float distNext = fabs(stats_shading_n[3] - stats_shading_s[3]);
            bool inf = (stats_shading_s[3] > stats_shading_n[3]);
			float expRes;
            if (distCur < distNext) {
                expRes = expCur;
            } else {
                expRes = expNext;
            }
            expArg = expIt;
            float diffAnglesSh = fabs(angle_saturation_s - angle_saturation_n);
            if (diffAnglesSh < minDiffAnglesSh) {
                minDiffAnglesSh = diffAnglesSh;
                argminDiffAnglesSh = absIt;
            }
            multiply(newAlbedo, newShading, newImage);
            string imageFilename (prefix), shadingFilename(prefix);
			imageFilename += "_res_"; shadingFilename += "_shading_";
			imageFilename += "_"; shadingFilename += "_";
			imageFilename += to_string(absIt); shadingFilename += to_string(absIt);
			imageFilename += "_"; shadingFilename += "_";
			imageFilename += to_string(int(100 * slope)); shadingFilename += to_string(int(100 * slope));
			imageFilename += "_"; shadingFilename += "_";
			imageFilename += to_string(int(100 * abs)); shadingFilename += to_string(int(100 * abs));
			imageFilename += "_"; shadingFilename += "_";
			imageFilename += to_string(int(100 * expRes)); shadingFilename += to_string(int(100 * expRes));
			imageFilename += ".png"; shadingFilename += ".png";
            _save(imageFilename, newImage);
			_save(shadingFilename, newShading);
            fprintf(logFile, "It. %2i (break = %6.3f) : Minimum found in exp = %5.1f : ", absIt, abs, expRes);
            fprintf(logFile, "mean = %9.5f, sd = %9.5f, skew = %9.5f, angle = %9.5f",
                    stats_shading_n[1], stats_shading_n[2], stats_shading_n[3], angle_saturation_n);
            if (!inf || !sup) {
                fprintf(logFile, " *");
            }
            fprintf(logFile, "\n");
        }
        fprintf(logFile, "\n");
        fprintf(logFile, "Optimal suggested result using shading saturation : %i\n", argminDiffAnglesSh);
        fclose(logFile);
    } else {
        std::cout << "generateOptimalShadingsNeg() : error : couldn't open a log file" << std::endl;
    }
}


void Method::multiply (Mat & albedo, Mat & shading, Mat & result) {
	for (int y = 0 ; y < albedo.rows ; y++) {
		for (int x = 0 ; x < albedo.cols ; x++) {
			result.at<Vec4f>(y, x)[0] = clamp(0, albedo.at<Vec4f>(y, x)[0] * shading.at<Vec4f>(y, x)[0], 1);
			result.at<Vec4f>(y, x)[1] = clamp(0, albedo.at<Vec4f>(y, x)[1] * shading.at<Vec4f>(y, x)[1], 1);
			result.at<Vec4f>(y, x)[2] = clamp(0, albedo.at<Vec4f>(y, x)[2] * shading.at<Vec4f>(y, x)[2], 1);
			result.at<Vec4f>(y, x)[3] = albedo.at<Vec4f>(y, x)[3] * shading.at<Vec4f>(y, x)[3];
		}
	}
}

