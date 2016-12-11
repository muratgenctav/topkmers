/*
 * topkmers.cpp
 *
 * Program main
 *
 *  Created on: Dec 6, 2016
 *      Author: muratgenctav
 */

#include <cstdio>
#include <vector>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "TopKmers.h"

using namespace std;

#define MAX_K 30						// max value for k
#define MAX_TOPCOUNT 25					// max value for topcount
#define MAX_NTHREADS 4					// max value for the number of threads
#define MAX_MAP_SIZE 10000000			// max size for the hashmap
#define OPT_HELP "--help"				// option help
#define OPT_INPUT "--input"				// option input
#define OPT_KMERLENGTH "--kmerlength"	// option k-mer length
#define OPT_TOPCOUNT "--topcount"		// option topcount
#define OPT_NTHREADS "--numthreads"		// option number of threads

/**
 * Function to show help
 */
void showHelp(string prog) {
	cerr << "Usage: " << endl
		 << prog << " " << OPT_HELP << endl
		 << "  Show this help message." << endl << endl
		 << prog << " " << OPT_INPUT << " <PATH> " << OPT_KMERLENGTH << " <LENGTH> [--other_options]" << endl
		 << "  Find and list top k-mers that appear in the specified fastq sequence file." << endl
		 << "  Options:" << endl
		 << "  " << OPT_INPUT << " <PATH>\t\tPath to fastq file (required)" << endl
		 << "  " << OPT_KMERLENGTH << " <LENGTH>\t\tLength 'k' of kmers up to " << MAX_K << " (required)" << endl
		 << "  " << OPT_TOPCOUNT << " <COUNT>\t\tNumber of top kmers to be listed up to " << MAX_TOPCOUNT << " (defaults to 1)" << endl
		 << "  " << OPT_NTHREADS << " <NUMTHREADS>\tNumber of threads up to " << MAX_NTHREADS << " (defaults to 1)" << endl;
}

/**
 * Function to check ranges of input parameters
 *
 * @param-in k: length of k-mer
 * @param-in nTopKmers: number of top k-mers to be listed
 * @param-in nThreads: number of counting threads
 * @return : 1 when check fails, 0 when OK
 */
int constraintsCheck(int k, int nTopKmers, int nThreads) {
	if(k < 1 || k > MAX_K) {
		cerr << "Please specify a positive argument for option " << OPT_KMERLENGTH << " less than or equal to " << MAX_K << endl;
		return 1;
	}
	if(nThreads < 1 || nThreads > MAX_NTHREADS) {
		cerr << "Please specify a positive argument for option " << OPT_NTHREADS << " less than or equal to " << MAX_NTHREADS << endl;
		return 1;
	}
	if(nTopKmers < 1 || nTopKmers > MAX_TOPCOUNT) {
		cerr << "Please specify a positive argument for option " << OPT_TOPCOUNT << " less than or equal to " << MAX_TOPCOUNT << endl;
		return 1;
	}
	return 0;
}

/**
 * Function to ensure consistency of input parameters
 *
 * @param-in k: length of k-mer
 * @param-in/out nTopKmers: number of top k-mers to be listed
 * @param-in/out nThreads: number of counting threads
 */
void consistencyCheck(int k, int &nTopKmers, int &nThreads) {
	if(log2(nTopKmers)/2 > (double) k) {
		int critical = (int)  pow(4,k);
		cout << "Warning: Number of top k-mers cannot exceed " << critical << " for the specified value of k = " << k << endl;
		nTopKmers = critical;
		cout << "Listing top " << nTopKmers << " k-mers instead." << endl;
	}
	if(log2(nThreads)/2 > (double) k) {
		int critical = (int)  pow(4,k);
		cout << "Warning: Number of threads cannot exceed " << critical << " for the specified value of k = " << k << endl;
		nThreads = critical;
		cout << "Using " << nThreads << " threads instead." << endl;
	}
}

/**
 * Program main function
 *
 * Processes input arguments, creates an instance of TopKmers class,
 * calls its method to compute top-kmers and prints returned list
 */
int main(int argc, char *argv[]) {
	// number of arguments check
    if(argc < 2) {
        showHelp(argv[0]);
        return 1;
    }
    // default values
    string inFile = "";
    int k = 0;
    int nTopKmers = 1;
    int nThreads = 1;
    // process arguments
    for(int i = 1; i < argc; i++) {
    	string arg = argv[i];
    	if(arg == OPT_HELP) {
    		showHelp(argv[0]);
    		return 0;
    	} else if(arg == OPT_INPUT) {
    		if(i+1 < argc) {
    			inFile = argv[++i];
    		} else {
    			cerr << OPT_INPUT << " option requires a path argument." << endl;
    			return 1;
    		}
    	} else if(arg == OPT_KMERLENGTH) {
    		if(i+1 < argc) {
    			if(!(stringstream(argv[++i]) >> k)) {
    				cerr << "Invalid argument for option " << OPT_KMERLENGTH << endl;
    				showHelp(argv[0]);
    				return 1;
    			}
    		} else {
    			cerr << OPT_KMERLENGTH << " option requires a number argument." << endl;
    			return 1;
    		}
    	} else if(arg == OPT_TOPCOUNT) {
    		if(i+1 < argc) {
    			if(!(stringstream(argv[++i]) >> nTopKmers)) {
    				cerr << "Invalid argument for option " << OPT_TOPCOUNT << endl;
    				showHelp(argv[0]);
    				return 1;
    			}
    		} else {
    			cerr << OPT_TOPCOUNT << " option requires a number argument." << endl;
    			return 1;
    		}
    	} else if(arg == OPT_NTHREADS) {
    		if(i+1 < argc) {
    			if(!(stringstream(argv[++i]) >> nThreads)) {
    				cerr << "Invalid argument for option " << OPT_NTHREADS << endl;
    				showHelp(argv[0]);
    				return 1;
    			}
    		} else {
    			cerr << OPT_NTHREADS << " option requires a number argument." << endl;
    			return 1;
    		}
    	} else {
    		cerr << "Invalid option name: " << arg << endl;
    		showHelp(argv[0]);
    		return 1;
    	}
    }
    // constraints check
    if(constraintsCheck(k,nTopKmers,nThreads)) {
    	return 1;
    }
    // consistency check
    consistencyCheck(k,nTopKmers,nThreads);
    // input file status check
    if(inFile == "") {
    	cerr << "Please specify a valid input file." << endl;
    	showHelp(argv[0]);
    	return 1;
    }
    ifstream infile(inFile);
    if(!infile.good()) {
        cerr << "Cannot open file \"" << inFile << "\"" << endl;
        infile.close();
        return 1;
    }
    infile.close();

    // count k-mers
    time_t start,finish;
    time(&start);
    // -
    TopKmers tk(inFile,k,nTopKmers,nThreads,MAX_MAP_SIZE);
    const vector<pair<string,unsigned int>> & topKmers = tk.getTopKmers();
    // -
    time(&finish);
    printf ("Completed in %.2lf seconds.\n", difftime(finish,start));


    // print top k-mers
    for(auto it = topKmers.cbegin(); it != topKmers.cend(); it++) {
        cout << it->first << " : " << it->second << endl;
    }

    return 0;
}
