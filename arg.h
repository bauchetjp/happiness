#ifndef _ARG_H_
#define _ARG_H_

#include <string>

using namespace std;

enum ColorDistance {
	RGB_DISTANCE,
	LAB_DISTANCE
};


enum MethodType {
	QUILTING
};


typedef struct {
	enum ColorDistance distance;
	enum MethodType type;
	union {
		struct {
			int quiltingIterations;
			int quiltingRadius;
		};
	};
	int trees;
	int recursions;
	bool automaticShading;
	float alpha;
	float beta;
	float gamma;
	string prefix;
} MethodParams;


typedef struct {
	string sourceBasename;
	string sourceNormalMapFilename;
	string sourceAlbedoFilename;
	string sourceShadingFilename;
	string targetBasename;
	string targetNormalMapFilename;
	string targetShadingFilename;
 	enum ColorDistance distance;
	MethodParams mparams;
} GeneralParams;


void setDefaultParameters (GeneralParams & gparams);

void usage ();

bool interpreter (int argc, char *argv[], GeneralParams & gparams);

void printParameters (GeneralParams & gparams);

#endif
