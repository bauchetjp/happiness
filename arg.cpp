#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "arg.h"


using namespace std;

void setDefaultParameters (ReconstructionParams & params) {
	params.sourceName = "";
	params.targetName = "";
	params.iterations = 3;
	params.radius = 10;
	params.distance = RGB_DISTANCE;
	params.automaticShading = true;
	params.alpha = 1.0f;
	params.beta = 0.0f;
	params.gamma = 1.0f;
	params.trees = 16;
	params.recursions = 1024;
}


void usage () {
	cout << "How to use this program : ./main [options] " << endl;
	cout << "where available options are : " << endl;
	cout << "    --source [object name] [AMBIENT | WEST | EAST | FRONT]" << endl;
	cout << "    --target [object name] [AMBIENT | WEST | EAST | FRONT]" << endl;
	cout << "    --iterations n" << endl;
	cout << "    --radius r" << endl;
	cout << "    --distance [RGB | LAB]" << endl;
	cout << "    --kdtrees-parameters trees recursions" << endl;
	cout << "    --shading-parameters [AUTO | MANUAL [alpha beta gamma]]" << endl;
	cout << endl;
	cout << "Please report bugs to bauchet.jeanphilippe@gmail.com." << endl;
}


bool interpreter (int argc, char *argv[], ReconstructionParams & params) {
	setDefaultParameters (params);
	int r = 1;
	while (r < argc) {
		if (!strcmp(argv[r], "--source") && (r + 2 < argc)) {
			params.sourceName = "./images/";
			params.sourceName += argv[r + 1];
			params.sourceName += "_";
			if (!strcmp(argv[r + 2], "AMBIENT")) {
				params.sourceName += "ambient";
			} else if (!strcmp(argv[r + 2], "WEST")) {
				params.sourceName += "west";
			} else if (!strcmp(argv[r + 2], "EAST")) {
				params.sourceName += "east";
			} else if (!strcmp(argv[r + 2], "FRONT")) {
				params.sourceName += "front";
			} else {
				usage();
				return false;
			}
			r += 3;
		} else if (!strcmp(argv[r], "--target") && (r + 2 < argc)) {
			params.targetName = "./images/";
			params.targetName += argv[r + 1];
			params.targetName += "_";
			if (!strcmp(argv[r + 2], "AMBIENT")) {
				params.targetName += "ambient";
			} else if (!strcmp(argv[r + 2], "WEST")) {
				params.targetName += "west";
			} else if (!strcmp(argv[r + 2], "EAST")) {
				params.targetName += "east";
			} else if (!strcmp(argv[r + 2], "FRONT")) {
				params.targetName += "front";
			} else {
				usage();
				return false;
			}
			r += 3;
		} else if (!strcmp(argv[r], "--iterations") && (r + 1 < argc)) {
			params.iterations = atoi(argv[r + 1]);
			r += 2;
		} else if (!strcmp(argv[r], "--radius") && (r + 1 < argc)) {
			params.radius = atoi(argv[r + 1]);
			r += 2;
		} else if (!strcmp(argv[r], "--distance") && (r + 1 < argc)) {
			if (!strcmp(argv[r + 1], "RGB")) {
				params.distance = RGB_DISTANCE;
			} else if (!strcmp(argv[r + 1], "LAB")) {
				params.distance = LAB_DISTANCE;
			} else {
				usage();
				return false;
			}
			r += 2;
		} else if (!strcmp(argv[r], "--kdtrees-parameters") && (r + 2 < argc)) {
			params.trees = atoi(argv[r + 1]);
			params.recursions = atoi(argv[r + 2]);
			r += 3;
		} else if (!strcmp(argv[r], "--shading-parameters")) {
			if (!strcmp(argv[r + 1], "AUTO") && r + 1 < argc) {
				params.automaticShading = true;
				r += 2;
			} else if (!strcmp(argv[r + 1], "MANUAL") && r + 4 < argc) {
				params.automaticShading = false;
				params.alpha = atof(argv[r + 2]);
				params.beta  = atof(argv[r + 3]);
				params.gamma = atof(argv[r + 4]);
				r += 5;
			} else {
				usage();
				return false;
			}
		} else if (!strcmp(argv[r], "--help")) {
			usage();
			return false;
		} else {
			usage();
			return false;
		}
	}
	return true;
}


void printParameters (ReconstructionParams & params) {
	cout << "Source : " << params.sourceName << endl;
	cout << "Target : " << params.targetName << endl;
	cout << "Iterations : " << params.iterations << endl;
	cout << "Radius : " << params.radius << endl;
	if (params.distance == RGB_DISTANCE) {
		cout << "Distance : RGB" << endl;
	} else {
		cout << "Distance : LAB" << endl;
	}
	if (params.automaticShading) {
		cout << "Automatic shading : true" << endl;
	} else {
		cout << "Automatic shading : false" << endl;
	}
	cout << "Alpha : " << params.alpha << endl;
	cout << "Beta : " << params.beta << endl;
	cout << "Gamma : " << params.gamma << endl;
	cout << "Trees : " << params.trees << endl;
	cout << "Recursions : " << params.recursions << endl;
}
