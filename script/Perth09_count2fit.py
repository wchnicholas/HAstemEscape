#!/usr/bin/python
import os
import sys
import glob
import string
import numpy as np
from Bio import SeqIO
from collections import Counter, defaultdict

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def CsvWithHeader2Hash(fitfile):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  header = []
  for line in infile.xreadlines():
    countline += 1
    line = line.rstrip().rsplit(",")
    if countline == 1: header = line; continue
    mut = line[0]
    H[mut] = {}
    for i in range(1,len(line)): H[mut][header[i]] = line[i]
  infile.close()
  return H

def nuc_to_aa_dict(nuc_count_dict):
  aa_count_dict = {}
  for pos in nuc_count_dict.keys():
    pos_count_dict = defaultdict(int)
    wtaa = translation(nuc_count_dict[pos]['wildtype'])
    ID = wtaa+pos.replace('(HA2)','')+'-'+'(HA2)' if '(HA2)' in pos else wtaa+pos
    for codon in nuc_count_dict[pos].keys():
      if codon=='wildtype': continue
      aa = translation(codon)
      pos_count_dict[aa] += int(nuc_count_dict[pos][codon])+1
      pos_count_dict['total'] += int(nuc_count_dict[pos][codon])+1
      if aa==wtaa:
        pos_count_dict['wtaa'] += int(nuc_count_dict[pos][codon])+1
    aa_count_dict[ID] = pos_count_dict
  return aa_count_dict

def count_to_freq(count_dict, sampleID, residue, aa):
  return float(count_dict[sampleID][residue][aa])/float(count_dict[sampleID][residue]['total'])

def compile_out_Perth09(count_dict, outfile):
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_']
  residues = list(set([residue for ID in count_dict.keys() for residue in count_dict[ID].keys()]))
  print "writiug %s" % outfile
  outfile = open(outfile,'w')
  outfile.write("\t".join(map(str,['mutID', 'fit_L1_mock', 'fit_L2_mock', 'fit_L3_mock',
                                   'fit_L2_FI6v3_5ug', 'fit_L2_FI6v3_12ug', 'fit_L3_FI6v3_13ug',
                                   'fit_L1_FI6v3_15ug', 'fit_L3_FI6v3_15ug',
                                   'fit_mock', 'fit_FI6v3_15ug', 'esp_FI6v3_15ug']))+"\n")
  for residue in residues:
    for aa in aas:
      mutID = residue.rsplit('-')[0]+aa+'-(HA2)' if '(HA2)' in residue else residue+aa
      total_mutDNA1 = count_dict['mutDNA-1'][residue]['total']
      total_mutDNA2 = count_dict['mutDNA-2'][residue]['total']
      total_mutDNA3 = count_dict['mutDNA-3'][residue]['total']

      #DEFINTE WT FITNESS
      wt_count_mutDNA1 = count_dict['mutDNA-1'][residue]['wtaa']
      wt_count_mutDNA2 = count_dict['mutDNA-2'][residue]['wtaa']
      wt_count_mutDNA3 = count_dict['mutDNA-3'][residue]['wtaa']
      wt_freq_mutDNA = float(wt_count_mutDNA1+wt_count_mutDNA2+wt_count_mutDNA3)/float(total_mutDNA1+total_mutDNA2+total_mutDNA3)

      wt_fit_L1_mock = count_to_freq(count_dict, 'L1mock', residue, 'wtaa')/wt_freq_mutDNA
      wt_fit_L2_mock = count_to_freq(count_dict, 'L2mock', residue, 'wtaa')/wt_freq_mutDNA
      wt_fit_L3_mock = count_to_freq(count_dict, 'L3mock', residue, 'wtaa')/wt_freq_mutDNA
      wt_fit_L2_FI6v3_5ug = count_to_freq(count_dict,'L2-FI6v3-5ug-ml', residue, 'wtaa')/wt_freq_mutDNA
      wt_fit_L2_FI6v3_12ug = count_to_freq(count_dict,'L2-FI6v3-12ug-ml', residue, 'wtaa')/wt_freq_mutDNA
      wt_fit_L3_FI6v3_13ug = count_to_freq(count_dict,'L3-FI6v3-13ug-ml', residue, 'wtaa')/wt_freq_mutDNA
      wt_fit_L1_FI6v3_15ug = count_to_freq(count_dict,'L1-FI6v3-15ug-ml', residue, 'wtaa')/wt_freq_mutDNA
      wt_fit_L3_FI6v3_15ug = count_to_freq(count_dict,'L3-FI6v3-15ug-ml', residue, 'wtaa')/wt_freq_mutDNA
      
      #COMPUTE MUTANT FITNESS
      count_mutDNA1 = count_dict['mutDNA-1'][residue][aa]
      count_mutDNA2 = count_dict['mutDNA-2'][residue][aa]
      count_mutDNA3 = count_dict['mutDNA-3'][residue][aa]
      freq_mutDNA = float(count_mutDNA1+count_mutDNA2+count_mutDNA3)/float(total_mutDNA1+total_mutDNA2+total_mutDNA3)

      fit_L1_mock = count_to_freq(count_dict, 'L1mock', residue, aa)/freq_mutDNA/wt_fit_L1_mock
      fit_L2_mock = count_to_freq(count_dict, 'L2mock', residue, aa)/freq_mutDNA/wt_fit_L2_mock
      fit_L3_mock = count_to_freq(count_dict, 'L3mock', residue, aa)/freq_mutDNA/wt_fit_L3_mock
      fit_L2_FI6v3_5ug = count_to_freq(count_dict,'L2-FI6v3-5ug-ml', residue, aa)/freq_mutDNA/wt_fit_L2_FI6v3_5ug
      fit_L2_FI6v3_12ug = count_to_freq(count_dict,'L2-FI6v3-12ug-ml', residue, aa)/freq_mutDNA/wt_fit_L2_FI6v3_12ug
      fit_L3_FI6v3_13ug = count_to_freq(count_dict,'L3-FI6v3-13ug-ml', residue, aa)/freq_mutDNA/wt_fit_L3_FI6v3_13ug
      fit_L1_FI6v3_15ug = count_to_freq(count_dict,'L1-FI6v3-15ug-ml', residue, aa)/freq_mutDNA/wt_fit_L1_FI6v3_15ug
      fit_L3_FI6v3_15ug = count_to_freq(count_dict,'L3-FI6v3-15ug-ml', residue, aa)/freq_mutDNA/ wt_fit_L3_FI6v3_15ug

      fit_mock = np.mean([fit_L1_mock, fit_L2_mock, fit_L3_mock])
      fit_FI6v3_15ug = np.mean([fit_L1_FI6v3_15ug, fit_L3_FI6v3_15ug])
      fitcutoff = 0.5
      esp_FI6v3_15ug = fit_FI6v3_15ug/fit_mock if fit_mock >= fitcutoff else -1
      outfile.write("\t".join(map(str,[mutID, fit_L1_mock, fit_L2_mock, fit_L3_mock,
                                       fit_L2_FI6v3_5ug, fit_L2_FI6v3_12ug, fit_L3_FI6v3_13ug, 
                                       fit_L1_FI6v3_15ug, fit_L3_FI6v3_15ug,
                                       fit_mock, fit_FI6v3_15ug, esp_FI6v3_15ug]))+"\n")
  outfile.close()

def main():
  filenames = glob.glob('Bloom_data/Perth09/*.csv')
  outfile   = 'result/Fitness_Perth09.tsv'
  count_dict = {}
  for filename in filenames:
    ID = filename.rsplit('/')[-1].rsplit('_')[0] 
    nuc_count_dict = CsvWithHeader2Hash(filename)
    aa_count_dict  = nuc_to_aa_dict(nuc_count_dict)
    count_dict[ID] = aa_count_dict
  compile_out_Perth09(count_dict, outfile)

if __name__ == "__main__":
  main()
