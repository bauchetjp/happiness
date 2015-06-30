#ifndef _ARG_H_
#define _ARG_H_

#include <string>

using namespace std;

enum ColorDistance {
	RGB_DISTANCE,
	LAB_DISTANCE
} ;

typedef struct {
	string sourceName;
	string targetName;
	int iterations;	
	int radius;
 	enum ColorDistance distance;
	bool automaticShading;
	float alpha;
	float beta;
	float gamma;
	int trees;
	int recursions;
} ReconstructionParams;


void setDefaultParameters (ReconstructionParams & params);

void usage ();

bool interpreter (int argc, char *argv[], ReconstructionParams & params);

void printParameters (ReconstructionParams & params);

#endif
