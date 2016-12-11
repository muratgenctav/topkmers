/*
 * TopKmers.h
 *
 * Class for finding most frequent k-mers
 *
 *  Created on: Dec 5, 2016
 *      Author: muratgenctav
 */

#ifndef TOPKMERS_H_
#define TOPKMERS_H_

#include <string>
#include <vector>
#include <unordered_map>

typedef unsigned long long int kmer_key_t;

using namespace std;

class TopKmers {
public:
	/**
	 * Constructor method
	 *
	 * @param-in fileName: the path of the FASTQ file
	 * @param-in k: the length of k-mers to be counted
	 * @param-in nTopKmers: number of top k-mers to be listed
	 * @param-in nThreads: number of threads, each counting a distinct partition of the all possible k-mers
	 * @param-in maxMapSize: constraint for the maximum allowed size of the hashmap(s) to use memory efficiently
	 */
	TopKmers(string fileName, int k, int nTopKmers, int nThreads, unsigned int maxMapSize = 1000000);

	/**
	 * Public method that extracts (only at first call) and returns top k-mers as (k-mer,frequency) pairs.
	 * In multi-threaded mode, partitions the k-mer space, creates counter threads and merges the results.
	 *
	 * @return const reference to the member data topKmers which stores top k-mers as (k-mer, frequency) pairs.
	 */
	const vector<pair<string,unsigned int>> & getTopKmers();

	// Destructor
	virtual ~TopKmers();

private:
	/**
	 * Thread that counts a partition of the k-mer space and returns top k-mers of the sub-space.
	 * In single thread mode, explores the entire space of possible k-mers.
	 *
	 * @param-out mostFrequentKmers: Output parameter to return the list of top k-mers
	 * @param-in threadID: Counter thread identifier
	 */
	void computeTopKmers(vector<pair<string,unsigned int>> &mostFrequentKmers, const int threadID = 0);

	/**
	 * Helper method to process a sequence
	 *
	 * @param-in sequence: DNA sequence to process
	 * @param-in/out countMap: count-map to be updated
	 * @param-in startIdx: starting index of the k-mer sub-space
	 * @param-in endIdx: ending index of the k-mer sub-space
	 */
	void processSeq(const string &sequence, unordered_map<kmer_key_t,unsigned int> &countMap, const kmer_key_t startIdx, const kmer_key_t endIdx);

	/**
	 * Helper method to extract most frequent k-mers from a count-map
	 *
	 * @param-out mostFrequentKmers: Output parameter to return the list of top k-mers
	 * @param-in countMap: count-map that maps k-mers to their counts
	 */
	void mostFreqKmers(vector<pair<string,unsigned int>> &mostFrequentKmers,const unordered_map<kmer_key_t,unsigned int> &countMap);

	/**
	 * Helper method to merge results from multiple threads (called when multi-thread mode)
	 * Merges them into the data member topKmers
	 *
	 * @param-in results: array of top k-mer lists output by threads
	 */
	void mergeMultipleResults(vector<vector<pair<string,unsigned int>>> const &results);

	// Class data members
	const string fileName;		// Path of the file to be processed
	const int k;				// The length of k-mers to be counted
	const int nTopKmers;		// The number of top k-mers expected to be listed
	const int nThreads;			// Number of threads that counts a partition of the all possible k-mers
	unsigned int maxMapSize;	// Constraint for the maximum size of the hashmap(s) allowed to use memory efficiently
	vector<pair<string,unsigned int>> topKmers; // List of extracted k-mers stored as (k-mer,frequency) pairs
};

#endif /* TOPKMERS_H_ */
