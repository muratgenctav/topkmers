/*
 * SeqFileScanner.h
 *
 * FASTQ scanner class
 *
 *  Created on: Dec 6, 2016
 *      Author: muratgenctav
 */

#ifndef SEQFILESCANNER_H_
#define SEQFILESCANNER_H_

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class SeqFileScanner {

private:
    string fileName;
    ifstream seqFile;

public:
    /**
     * Constructor
     *
     * @param-in fName: the path of the FASTQ file
     */
    SeqFileScanner(string fName)
    : fileName(fName) {}

    /**
     * Method for reading the next sequence
     *
     * @param-out seqBuff: buffer string to read the next sequence into
     * @return : 1 when read is successful, 0 otherwise
     */
    int readNextSequence(string &seqBuff) {
        if(seqFile.ignore(numeric_limits<streamsize>::max(),'\n')){
            if(getline(seqFile,seqBuff)) {
                seqFile.ignore(numeric_limits<streamsize>::max(),'\n');
                seqFile.ignore(numeric_limits<streamsize>::max(),'\n');
                return 1;
            }
        }
        return 0;
    }

    /**
     * Method for opening FASTQ file
     *
     * @return : 1 when the file is successfully opened, 0 otherwise
     */
    int open() {
        seqFile.open(fileName);
        if(seqFile.good()) {
            return 1;
        }
        return 0;
    }

    /**
     * Method for closing FASTQ file
     */
    void close() {
        if(seqFile.is_open()) {
            seqFile.close();
        }
    }

    /**
     * Destructor
     */
    ~SeqFileScanner() {
        close();
    }
};

#endif /* SEQFILESCANNER_H_ */
