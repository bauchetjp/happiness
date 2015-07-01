#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "arg.h"


using namespace std;

void setDefaultParameters (GeneralParams & gparams) {
	gparams.sourceBasename = "";
	gparams.sourceNormalMapFilename = "";
	gparams.sourceAlbedoFilename = "";
	gparams.sourceShadingFilename = "";
	gparams.targetBasename = "";
	gparams.targetNormalMapFilename = "";
	gparams.targetShadingFilename = "";
	gparams.mparams.distance = RGB_DISTANCE;
	gparams.mparams.type = QUILTING;
	gparams.mparams.quiltingIterations = 3;
	gparams.mparams.quiltingRadius = 10;
	gparams.mparams.automaticShading = true;
	gparams.mparams.alpha = 1.0f;
	gparams.mparams.beta = 0.0f;
	gparams.mparams.gamma = 1.0f;
	gparams.mparams.trees = 16;
	gparams.mparams.recursions = 1024;
}


void usage () {
	cout << "How to use this program : ./main [options] " << endl;
	cout << "where available options are : " << endl;
	cout << "    --source [object name] [AMBIENT | WEST | EAST | FRONT]" << endl;
	cout << "    --target [object name] [AMBIENT | WEST | EAST | FRONT]" << endl;
	cout << "    --method [QUILTING [iterations radius]]" << endl;
	cout << "    --distance [RGB | LAB]" << endl;
	cout << "    --kdtrees-parameters trees recursions" << endl;
	cout << "    --shading-parameters [AUTO | MANUAL [alpha beta gamma]]" << endl;
	cout << endl;
	cout << "Please report bugs to bauchet.jeanphilippe@gmail.com." << endl;
}


bool interpreter (int argc, char *argv[], GeneralParams & params) {
	setDefaultParameters (params);
	int r = 1;
	while (r < argc) {
		if (!strcmp(argv[r], "--source") && (r + 2 < argc)) {
			params.sourceBasename = "./images/";
			params.sourceBasename += argv[r + 1];
			params.sourceBasename += "_";
			if (!strcmp(argv[r + 2], "AMBIENT")) {
				params.sourceBasename += "ambient";
			} else if (!strcmp(argv[r + 2], "WEST")) {
				params.sourceBasename += "west";
			} else if (!strcmp(argv[r + 2], "EAST")) {
				params.sourceBasename += "east";
			} else if (!strcmp(argv[r + 2], "FRONT")) {
				params.sourceBasename += "front";
			} else {
				usage();
				return false;
			}
			params.sourceNormalMapFilename = params.sourceBasename;
			params.sourceAlbedoFilename = params.sourceBasename;
			params.sourceShadingFilename = params.sourceBasename;
			params.sourceNormalMapFilename += "_n.txt";
			params.sourceAlbedoFilename += "_a.txt";
			params.sourceShadingFilename += "_s.txt";
			r += 3;
		} else if (!strcmp(argv[r], "--target") && (r + 2 < argc)) {
			params.targetBasename = "./images/";
			params.targetBasename += argv[r + 1];
			params.targetBasename += "_";
			if (!strcmp(argv[r + 2], "AMBIENT")) {
				params.targetBasename += "ambient";
			} else if (!strcmp(argv[r + 2], "WEST")) {
				params.targetBasename += "west";
			} else if (!strcmp(argv[r + 2], "EAST")) {
				params.targetBasename += "east";
			} else if (!strcmp(argv[r + 2], "FRONT")) {
				params.targetBasename += "front";
			} else {
				usage();
				return false;
			}
			params.targetNormalMapFilename = params.targetBasename;
			params.targetShadingFilename = params.targetBasename;
			params.targetNormalMapFilename += "_n.txt";
			params.targetShadingFilename += "_s.txt";
			r += 3;
		} else if (!strcmp(argv[r], "--method") && (r + 1 < argc)) {
			if (!strcmp(argv[r + 1], "QUILTING") && (r + 3 < argc)) {
				params.mparams.type = QUILTING;
				params.mparams.quiltingIterations = atoi(argv[r + 2]);
				params.mparams.quiltingRadius = atoi(argv[r + 3]);
				r += 4;
			} else {
				usage();
				return false;
			}
		} else if (!strcmp(argv[r], "--distance") && (r + 1 < argc)) {
			if (!strcmp(argv[r + 1], "RGB")) {
				params.mparams.distance = RGB_DISTANCE;
			} else if (!strcmp(argv[r + 1], "LAB")) {
				params.mparams.distance = LAB_DISTANCE;
			} else {
				usage();
				return false;
			}
			r += 2;
		} else if (!strcmp(argv[r], "--kdtrees-parameters") && (r + 2 < argc)) {
			params.mparams.trees = atoi(argv[r + 1]);
			params.mparams.recursions = atoi(argv[r + 2]);
			r += 3;
		} else if (!strcmp(argv[r], "--shading-parameters")) {
			if (!strcmp(argv[r + 1], "AUTO") && r + 1 < argc) {
				params.mparams.automaticShading = true;
				r += 2;
			} else if (!strcmp(argv[r + 1], "MANUAL") && r + 4 < argc) {
				params.mparams.automaticShading = false;
				params.mparams.alpha = atof(argv[r + 2]);
				params.mparams.beta  = atof(argv[r + 3]);
				params.mparams.gamma = atof(argv[r + 4]);
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


void printParameters (GeneralParams & gparams) {
	cout << "Source normal map : " << gparams.sourceNormalMapFilename << endl;
	cout << "Target normal map : " << gparams.targetNormalMapFilename << endl;
	cout << "Iterations : " << gparams.mparams.quiltingIterations << endl;
	cout << "Radius : " << gparams.mparams.quiltingRadius << endl;
	if (gparams.mparams.distance == RGB_DISTANCE) {
		cout << "Distance : RGB" << endl;
	} else {
		cout << "Distance : LAB" << endl;
	}
	if (gparams.mparams.automaticShading) {
		cout << "Automatic shading : true" << endl;
	} else {
		cout << "Automatic shading : false" << endl;
	}
	cout << "Alpha : " << gparams.mparams.alpha << endl;
	cout << "Beta : " << gparams.mparams.beta << endl;
	cout << "Gamma : " << gparams.mparams.gamma << endl;
	cout << "Trees : " << gparams.mparams.trees << endl;
	cout << "Recursions : " << gparams.mparams.recursions << endl;
}
