#!usr/bin/python
import os
import sys
import glob
from math import log
from collections import defaultdict

def TsvWithHeader2Hash(fitfile):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  header = []
  for line in infile.xreadlines():
    countline += 1
    line = line.rstrip().rsplit("\t")
    if countline == 1: header = line; continue
    mut = line[0]
    H[mut] = {}
    for i in range(1,len(line)): H[mut][header[i]] = line[i]
  infile.close()
  return H

def format_heatmap(D_fit_dict, outfile):
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  WTresi = ['Q42','I45','D46','Q47','I48','N49','L52','T111']
  print "Writing %s" % outfile
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['mut1','mut2','fit_no_antibody','resist_10ug_CR9114','resist_2500ng_FI6v3'])+"\n")
  for resi1 in WTresi:
    for resi2 in WTresi:
      pos1  = int(resi1[1::])
      pos2  = int(resi2[1::])
      wtaa1 = resi1[0]
      wtaa2 = resi2[0]
      if pos1 >= pos2: continue
      for aa1 in aas:
        if wtaa1 == aa1: continue
        for aa2 in aas:
          if wtaa2 == aa2: continue
          m1  = resi1+aa1
          m2  = resi2+aa2
          Dmut = m1+'/'+m2
          fit_no_antibody = -1
          resist_10ug_CR9114 = -1
          resist_2500ng_FI6v3 = -1
          if Dmut in D_fit_dict.keys():
            fit_no_antibody = D_fit_dict[Dmut]['fit_no_antibody']
            if float(fit_no_antibody) > 0.5: 
              resist_10ug_CR9114  = D_fit_dict[Dmut]['resist_10ug_CR9114']
              resist_2500ng_FI6v3 = D_fit_dict[Dmut]['resist_2500ng_FI6v3']
          outfile.write("\t".join(map(str,[m1, m2, fit_no_antibody, resist_10ug_CR9114, resist_2500ng_FI6v3]))+"\n")
  outfile.close()

def main():
  file_esp = "result/Resist_D.tsv"
  file_out = "result/Heatmap_D.tsv"
  D_fit_dict = TsvWithHeader2Hash(file_esp)
  format_heatmap(D_fit_dict, file_out)

if __name__ == "__main__":
  main()
