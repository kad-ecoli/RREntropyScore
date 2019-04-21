#!/usr/bin/env python

"""
The program rrces.py calculates
the residue-residue contact entropy score.
The lists of contacts are read from Andriy's 
program output.
"""

import sys
import os
import os.path
import re
from RREntropyScore import RREntropyScore
from RRJsonManager import RRJsonManager

# global parameters
# input directory
INDIR = "./input/"
# output directory
OUTDIR = "./output/"
# CASP targets directory containing files with sequences in fasta format
SEQ_DIR = "./seqs/"

# read the length of the whole target from sequence
def getSeqLen(target):
  seqfile = SEQ_DIR + target[0:5] + ".seq.txt"  
  try:
    fh = open(seqfile, "r")
    length = 0;
    for l in fh:
      l = l.rstrip()
      if (re.match("^[A-Z]+$", l)):
        length = length + len(l)
    fh.close()
    return length
  except IOError:
    return None

# get list of models from directory
def getModels(indir):
  models = []
  tmp = {}
  for f in os.listdir(indir):
    if (f[5:7] == "RR" and f[-4:] == 'json'):
      tmp[f[0:12]] = 1
  models = tmp.keys()
  models.sort()
  return models

# process target
def processTarget(target):
  # read seq length
  prot_len = getSeqLen(target)
  #print prot_len
  rrs = RREntropyScore() # envoke RREntropyScore class
  jm = RRJsonManager() # envoke RRJsonManager class
  # read models
  indir = INDIR + target + "_contList"
  models = getModels(indir)
  # loop over rr_ranges
  for _range in ('L', 'ML'):
    targetfile = indir + "/" + target + "_contacts" + _range + ".json"
    targetConts = jm.reformatTargetConts(jm.readConts(targetfile))
    outfile = OUTDIR + "/" + target + "_" + _range + "rrces"
    if (os.path.isfile(outfile)):
      continue
    out = open(outfile, "w")
    out.write("Model,10p0,L5p0,L2p0,FLp0,10p05,L5p05,L2p05,FLp05\n")
    # loop over models
    for m in models:
      out.write(m + ",")
      modelfile = indir + "/" + m + "_" + _range + "cont.json"
      rr_data = jm.readConts(modelfile)
      # loop over probabilities
      for _prob in [ 0.0, 0.5 ]:
        # loop over list sizes
        for _rr_len in [ "10", "L5", "L2", "FL" ]:
          rr_conts = jm.filterRRcontacts(rr_data, targetConts, prob=_prob, rr_len=_rr_len)
          rrces = rrs.calcScore(prot_len, rr_conts)
          out.write(str(rrces))
          if ( _rr_len == 'FL' and _prob == 0.5 ):
            out.write("\n")
          else:
            out.write(",")
    out.close()



if __name__ == "__main__":
  target = sys.argv[1]
  processTarget(target)
