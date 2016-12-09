/*
 * TopKmers.cpp
 *
 * 	Implementation of class methods to find most frequent k-mers
 *
 *  Created on: Dec 5, 2016
 *      Author: muratgenctav
 */

#include "TopKmers.h"
#include <cmath>
#include <queue>
#include <thread>
#include "SeqFileScanner.h"

// Constructor method
TopKmers::TopKmers(string _fileName, int _k, int _nTopKmers, int _nThreads, unsigned int _maxMapSize)
: fileName(_fileName), k(_k), nTopKmers(_nTopKmers), nThreads(_nThreads) {
	maxMapSize = _maxMapSize/_nThreads;
	nSignChars = ceil(log2(_nThreads)/2);
}

// Public method that extracts (only at first call) and returns top k-mers as (k-mer,frequency) pairs
const vector<pair<string,unsigned int>> & TopKmers::getTopKmers() {
	if(topKmers.empty()) {
		// Compute top k-mers first
		if(nThreads == 1) {
			// one thread
			computeTopKmers(topKmers);
		} else {
			// multi-thread
			vector<vector<pair<string,unsigned int>>> topKmersOfPartition(nThreads);
			thread th[nThreads];
			for(int i = 0; i < nThreads; i++) {
				th[i] = thread(&TopKmers::computeTopKmers,this,ref(topKmersOfPartition[i]),i);
			}
			for(int i = 0; i < nThreads; i++) {
				th[i].join();
			}

			mergeMultipleResults(topKmersOfPartition);
		}
	}
	return topKmers;
}

// Computation method (thread)
void TopKmers::computeTopKmers(vector<pair<string,unsigned int>> &mostFrequentKmers, const int threadID) {
    // count k-mers
    string seq;
    unordered_map<string,unsigned int> countmap;
    SeqFileScanner scanner(fileName);
    if(!scanner.open()) {
        cerr << "Error: Thread-" << threadID << " couldn't open file." << endl;
        return;
    }
    unsigned int startKmerIdx = 0;
    unsigned int endKmerIdx = pow(4,nSignChars);
    unsigned int rangePerThread = endKmerIdx/nThreads;
    if(nThreads > 1) {
        startKmerIdx = threadID*rangePerThread;
        if(threadID < nThreads-1) {
            endKmerIdx = startKmerIdx+rangePerThread;
        }
    }
    while(scanner.readNextSequence(seq)) {
        processSeq(seq,countmap,startKmerIdx,endKmerIdx);
    }
    scanner.close();
    // determine top k-mers
    mostFreqKmers(mostFrequentKmers,countmap);
}

// External utility function to determine k-mer's partition
unsigned int seq2idx(const string &s,int cc) {
    unsigned int idx = 0;
    unsigned int factor = 1;
    for(int i = 0; i < cc; i++) {
        unsigned int code;
        switch(s[i]) {
            case 'A': code = 0; break;
            case 'T': code = 3; break;
            case 'G': code = 2; break;
            case 'C': code = 1; break;
            default : code = 0;
        }
        idx += code*factor;
        factor *= 4;
    }
    return idx;
}

// Helper method to process a sequence
void TopKmers::processSeq(const string &seq, unordered_map<string,unsigned int> &count, const unsigned int startIdx, const unsigned int endIdx) {
    unsigned int endSeq = seq.length()-k;
	for(unsigned int i = 0; i <= endSeq; i++) {
        string kmer = seq.substr(i,k);
        if(nThreads > 1) {
			unsigned int idx = seq2idx(kmer,nSignChars);
			if(idx < startIdx || idx >= endIdx) {
				// k-mer is not within the partition of the thread
				continue;
			}
		}
        if(count.size() > maxMapSize) {
            if(count.find(kmer) == count.end()) {
                // k-mer shows up after how long!
            	// ignore it
                continue;
            }
        }
        count[kmer]++;
    }
}

// External utility function object to compare two k-mers w.r.t. their frequency values
struct CompareFreq {
    bool operator()(const pair<string,unsigned int> &n1,const pair<string,unsigned int> &n2) {
        return n1.second > n2.second;
    }
};

// Helper method to extract most frequent k-mers from a count-map
void TopKmers::mostFreqKmers(vector<pair<string,unsigned int>> &mfKmers,const unordered_map<string,unsigned int> &count) {
    priority_queue<pair<string,unsigned int>,vector<pair<string,unsigned int>>,CompareFreq> pq;
    int qSize = 0;
    for(auto it = count.cbegin(); it != count.cend(); it++) {
        if(qSize < nTopKmers) {
            // fill up heap
            pq.push(pair<string,unsigned int>(it->first,it->second));
            qSize++;
        } else {
            // replace smallest if new element is greater
            if(it->second > pq.top().second) {
                pq.pop();
                pq.push(pair<string,unsigned int>(it->first,it->second));
            }
        }
    }
    // now pq stores top k-mers and their counts
    mfKmers.resize(qSize);
    for(int i = qSize-1; i >= 0; i--) {
        mfKmers[i] = pq.top();
        pq.pop();
    }
}

// Helper method to merge results from multiple threads into topKmers data member (called when multi-thread mode)
void TopKmers::mergeMultipleResults(vector<vector<pair<string,unsigned int>>> const &results) {
    priority_queue<pair<string,unsigned int>,vector<pair<string,unsigned int>>,CompareFreq> pq;
    // fill up heap
    int qSize = 0;
    for(int i = 0; i < nThreads; i++) {
        for(auto it = results[i].cbegin(); it != results[i].cend(); it++) {
        	if(qSize < nTopKmers) {
        		// fill up heap
				pq.push(pair<string,unsigned int>(it->first,it->second));
				qSize++;
        	} else {
				// replace smallest if new element is greater
				if(it->second > pq.top().second) {
					pq.pop();
					pq.push(pair<string,unsigned int>(it->first,it->second));
				}
        	}
        }
    }
    // now pq stores top k-mers and their counts
    topKmers.resize(qSize);
    for(int i = qSize-1; i >= 0; i--) {
        topKmers[i] = pq.top();
        pq.pop();
    }
}

TopKmers::~TopKmers() {}
