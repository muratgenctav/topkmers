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
    unordered_map<kmer_key_t,unsigned int> countmap;
    SeqFileScanner scanner(fileName);
    if(!scanner.open()) {
        cerr << "Error: Thread-" << threadID << " couldn't open file." << endl;
        return;
    }
    kmer_key_t startKmerIdx = 0;
    kmer_key_t endKmerIdx = 1;
    for(int i = 0; i < k; i++) {
    	endKmerIdx *= 4;
    }
    kmer_key_t rangePerThread = endKmerIdx/(kmer_key_t) nThreads;
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

// Utility function to encode k-mer as an integer
kmer_key_t seq2key(const string &seq) {
	kmer_key_t idx;
	switch(seq[0]) {
	case 'A':
		idx = 0; break;
	case 'T':
		idx = 1; break;
	case 'G':
		idx = 2; break;
	case 'C':
		idx = 3; break;
	default:
		return numeric_limits<unsigned long long>::max(); // err idx
	}
	for(size_t i = 1; i < seq.length(); i++ ){
		idx *= 4;
		switch(seq[i]) {
		case 'A':
			idx += 0; break;
		case 'T':
			idx += 1; break;
		case 'G':
			idx += 2; break;
		case 'C':
			idx += 3; break;
		default:
			return numeric_limits<unsigned long long>::max(); // err idx
		}
	}
	return idx;
}
string key2seq(kmer_key_t idx, int length) {
	char seq[length+1];
	for(int i = length-1; i >= 0; i--) {
		kmer_key_t charCode = idx & 3;
		switch(charCode) {
		case 0:
			seq[i] = 'A'; break;
		case 1:
			seq[i] = 'T'; break;
		case 2:
			seq[i] = 'G'; break;
		case 3:
			seq[i] = 'C'; break;
		}
		idx /= 4;
	}
	seq[length] = '\0';
	return string(seq);
}

// Helper method to process a sequence
void TopKmers::processSeq(const string &seq, unordered_map<kmer_key_t,unsigned int> &count, const kmer_key_t startIdx, const kmer_key_t endIdx) {
    unsigned int endSeq = seq.length()-k;
	for(unsigned int i = 0; i <= endSeq; i++) {
        string kmer = seq.substr(i,k);
        kmer_key_t key = seq2key(kmer);
        if(nThreads > 1) {
			if(key < startIdx || key >= endIdx) {
				// k-mer is not within the partition of the thread
				continue;
			}
		}
        if(count.size() > maxMapSize) {
            if(count.find(key) == count.end()) {
                // k-mer shows up after how long!
            	// ignore it
                continue;
            }
        }
        count[key]++;
    }
}

// External utility function object to compare two k-mers w.r.t. their frequency values
template<class T> struct CompareFreq {
    bool operator()(const pair<T,unsigned int> &n1,const pair<T,unsigned int> &n2) {
        return n1.second > n2.second;
    }
};

// Helper method to extract most frequent k-mers from a count-map
void TopKmers::mostFreqKmers(vector<pair<string,unsigned int>> &mfKmers,const unordered_map<kmer_key_t,unsigned int> &count) {
    priority_queue<pair<kmer_key_t,unsigned int>,vector<pair<kmer_key_t,unsigned int>>,CompareFreq<kmer_key_t>> pq;
    int qSize = 0;
    for(auto it = count.cbegin(); it != count.cend(); it++) {
        if(qSize < nTopKmers) {
            // fill up heap
            pq.push(pair<kmer_key_t,unsigned int>(it->first,it->second));
            qSize++;
        } else {
            // replace smallest if new element is greater
            if(it->second > pq.top().second) {
                pq.pop();
                pq.push(pair<kmer_key_t,unsigned int>(it->first,it->second));
            }
        }
    }
    // now pq stores top k-mers and their counts
    mfKmers.resize(qSize);
    for(int i = qSize-1; i >= 0; i--) {
    	const pair<kmer_key_t,unsigned int> &top = pq.top();
        mfKmers[i] = pair<string,unsigned int>(key2seq(top.first,k),top.second);
        pq.pop();
    }
}

// Helper method to merge results from multiple threads into topKmers data member (called when multi-thread mode)
void TopKmers::mergeMultipleResults(vector<vector<pair<string,unsigned int>>> const &results) {
    priority_queue<pair<string,unsigned int>,vector<pair<string,unsigned int>>,CompareFreq<string>> pq;
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
