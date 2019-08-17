#!/usr/bin/python
import os
import sys
import glob
import string
from Bio import SeqIO
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

def compile_fitdict(fitdict_HK68, fitdict_WSN, fitdict_SI06, fitdict_Mich15, fitdict_Perth09, aas,residues, outfile):
  print "writing: %s" % outfile
  outfile = open(outfile,'w')
  outfile.write("\t".join(['pos', 'aa', 'HK68_fit_no_antibody', 'HK68_fit_2ug_CR9114', 'HK68_fit_10ug_CR9114',
                           'HK68_fit_300ng_FI6v3', 'HK68_fit_2500ng_FI6v3',
                           'WSN_fit_mock', 'WSN_fit_CR9114_50ng', 'WSN_fit_CR9114_70ng', 'WSN_fit_CR9114_100ng', 
                           'WSN_fit_FI6v3_100ng', 'WSN_fit_FI6v3_200ng',
                           'SI06_fit_no_antibody', 'SI06_fit_300ng_FI6v3', 
                           'Mich15_fit_no_antibody', 'Mich15_fit_300ng_FI6v3', 
                           'Perth09_fit_mock', 'Perth09_fit_FI6v3_15ug'])+"\n")
  for residue in residues:
    for aa in aas:
      HK68_fit_no_antibody            = '1.0'
      HK68_fit_2ug_CR9114            = '1.0'
      HK68_fit_10ug_CR9114            = '1.0'
      HK68_fit_300ng_FI6v3            = '1.0'
      HK68_fit_2500ng_FI6v3            = '1.0'
      for ID_WSN in fitdict_WSN.keys():
        if '(HA2)' in ID_WSN:
          ID_WSN_HA2 = ID_WSN.rsplit('-')[0][1::]
          if ID_WSN_HA2 == residue+aa:
            WSN_fit_mock         = fitdict_WSN[ID_WSN]['fit_mock']
            WSN_fit_CR9114_50ng  = fitdict_WSN[ID_WSN]['fit_CR9114_50ng']
            WSN_fit_CR9114_70ng  = fitdict_WSN[ID_WSN]['fit_CR9114_70ng']
            WSN_fit_CR9114_100ng = fitdict_WSN[ID_WSN]['fit_CR9114_100ng']
            WSN_fit_FI6v3_100ng = fitdict_WSN[ID_WSN]['fit_FI6v3_100ng']
            WSN_fit_FI6v3_200ng = fitdict_WSN[ID_WSN]['fit_FI6v3_200ng']
      for ID_SI06 in fitdict_SI06.keys():
        if fitdict_SI06[ID_SI06]['resi']+fitdict_SI06[ID_SI06]['aa'] == residue+aa:
          SI06_fit_no_antibody = fitdict_SI06[ID_SI06]['fit_no_antibody']
          SI06_fit_300ng_FI6v3 = fitdict_SI06[ID_SI06]['fit_300ng_FI6v3']
      for ID_Mich15 in fitdict_Mich15.keys():
        if fitdict_Mich15[ID_Mich15]['resi']+fitdict_Mich15[ID_Mich15]['aa'] == residue+aa:
          Mich15_fit_no_antibody = fitdict_Mich15[ID_Mich15]['fit_no_antibody']
          Mich15_fit_300ng_FI6v3 = fitdict_Mich15[ID_Mich15]['fit_300ng_FI6v3']
      for ID_Perth09 in fitdict_Perth09.keys():
        if '(HA2)' in ID_Perth09:
          ID_Perth09_HA2 = ID_Perth09.rsplit('-')[0][1::]
          if ID_Perth09_HA2 == residue+aa:
            Perth09_fit_mock       = fitdict_Perth09[ID_Perth09]['fit_mock']
            Perth09_fit_FI6v3_15ug = fitdict_Perth09[ID_Perth09]['fit_FI6v3_15ug']
      for ID_HK68 in fitdict_HK68.keys():
        if fitdict_HK68[ID_HK68]['resi']+fitdict_HK68[ID_HK68]['aa'] == residue+aa:
          HK68_fit_no_antibody = fitdict_HK68[ID_HK68]['fit_no_antibody']
          HK68_fit_2ug_CR9114 = fitdict_HK68[ID_HK68]['fit_2ug_CR9114']
          HK68_fit_10ug_CR9114 = fitdict_HK68[ID_HK68]['fit_10ug_CR9114']
          HK68_fit_300ng_FI6v3 = fitdict_HK68[ID_HK68]['fit_300ng_FI6v3']
          HK68_fit_2500ng_FI6v3 = fitdict_HK68[ID_HK68]['fit_2500ng_FI6v3']
      outfile.write("\t".join([residue, aa,
                               HK68_fit_no_antibody, HK68_fit_2ug_CR9114, HK68_fit_10ug_CR9114,
                               HK68_fit_300ng_FI6v3, HK68_fit_2500ng_FI6v3,
                               WSN_fit_mock, WSN_fit_CR9114_50ng, WSN_fit_CR9114_70ng,
                               WSN_fit_CR9114_100ng, WSN_fit_FI6v3_100ng, WSN_fit_FI6v3_200ng,
                               SI06_fit_no_antibody, SI06_fit_300ng_FI6v3,
                               Mich15_fit_no_antibody, Mich15_fit_300ng_FI6v3,
                               Perth09_fit_mock, Perth09_fit_FI6v3_15ug])+"\n")
  outfile.close()

def main():
  outfile = 'result/Fit_compare.tsv'
  fitfile_HK68    = 'result/Resist_S.tsv'
  fitfile_WSN     = 'result/Fitness_WSN.tsv'
  fitfile_SI06    = 'result/Resist_SI06.tsv'
  fitfile_Mich15  = 'result/Resist_Mich15.tsv'
  fitfile_Perth09 = 'result/Fitness_Perth09.tsv'
  aas      = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  residues = ['42','45','46','47','48','49','52','111']
  fitdict_HK68    = TsvWithHeader2Hash(fitfile_HK68)
  fitdict_WSN     = TsvWithHeader2Hash(fitfile_WSN)
  fitdict_SI06    = TsvWithHeader2Hash(fitfile_SI06)
  fitdict_Mich15  = TsvWithHeader2Hash(fitfile_Mich15)
  fitdict_Perth09 = TsvWithHeader2Hash(fitfile_Perth09)
  compile_fitdict(fitdict_HK68, fitdict_WSN, fitdict_SI06, fitdict_Mich15, fitdict_Perth09, aas, residues, outfile)

if __name__ == "__main__":
  main()
