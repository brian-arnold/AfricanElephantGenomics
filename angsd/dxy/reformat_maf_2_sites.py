#!/usr/bin/python -tt

import re
import sys
import os
import gzip

def main():

	scaffold = sys.argv[1]

	fn = scaffold + "_All.mafs.gz"
	fh = gzip.open(fn, 'rt', encoding='utf-8')
	outFile = scaffold + ".sites"
	oh = open(outFile, 'w')
	for line in fh:
		line = line.strip()
		info = line.split()
		if line.startswith("scaffold"):
			print(info[0], info[1], file=oh)
	oh.close()	

if __name__ == '__main__':
  main()
