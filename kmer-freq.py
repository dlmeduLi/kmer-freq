#!/usr/bin/env python

from __future__ import print_function
import re
from optparse import OptionParser
import sys
import os
import os.path

def AnalyseKmer(read, ksize, kmerRe, dictKmer):
	for k in range(ksize):
		kmers = kmerRe.findall(read[k:])
		i = 0
		for kmer in kmers:
			if(not kmer in dictKmer):
				dictKmer[kmer] = [0] * len(read)
			offset = k + i * ksize
			dictKmer[kmer][offset] += 1
			i += 1

def WriteKmer(dictKmer, outFile):
	for kmer in dictKmer:
		outFile.write("%s,%s\n" % (kmer, ','.join(map(str, dictKmer[kmer]))))

def main():

	# parse the command line options

	usage = 'usage: %prog [options] input.fq KSIZE'
	parser = OptionParser(usage=usage, version='%prog 0.1.0')

	(options, args) = parser.parse_args()
	if(len(args) != 2):
		parser.print_help()
		sys.exit(0)

	inputFileName = args[0]
	ksize = int(args[1])

	kmerRegExp = '(.{' + str(ksize) + '})' 
	kmerRe = re.compile(kmerRegExp)
	
	if(not os.path.exists(inputFileName)):
		print('error: fastq file "', inputFileName, '"', ' doest not exist.')
		sys.exit(-1)

	baseFileName = os.path.splitext(os.path.basename(inputFileName))[0]
	outFileName = baseFileName + '.' + str(ksize) + 'mer'
	try:
		outFile = open(outFileName, 'w')
	except IOError:
		print('error: Create output file failed!')
		sys.exit(-1)

	# kmer frequency table format dictKmer = {"CGTCG" : [0, 1, 2, ..., readlen - 1]}

	dictKmer = {}
	lineNum = 0
	with open(inputFileName, 'r') as fqFile :
		for line in fqFile:
			if((lineNum + 3) % 4 == 0):
				AnalyseKmer(line, ksize, kmerRe, dictKmer)
			lineNum += 1	
		WriteKmer(dictKmer, outFile)

	fqFile.close()

if __name__ == '__main__':
	main()