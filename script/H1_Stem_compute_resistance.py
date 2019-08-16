#!/usr/bin/python
import os
import sys
import glob

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

def fit2escape(fit_dict, outfile, inputcutoff, WTresi, aas, fitcutoff):
  print "writing %s" % outfile
  all_S_muts = [resi+aa for resi in WTresi for aa in aas]
  WT_ipt_count = ''
  outfile = open(outfile,'w')
  outfile.write("\t".join(['mut','resi','aa','Input_count','fit_no_antibody',
                           'fit_300ng_FI6v3','resist_300ng_FI6v3'])+"\n")
  for mut in fit_dict.keys():
    resi  = mut[1:-1]
    aa    = mut[-1] 
    Input_count = float(fit_dict[mut]['Input_count'])
    rep1_no_antibody_fit    = float(fit_dict[mut]['rep1_no_antibody_fit'])
    rep2_no_antibody_fit    = float(fit_dict[mut]['rep2_no_antibody_fit'])
    rep1_300ng_FI6v3_fit    = float(fit_dict[mut]['rep1_300ng_FI6v3_fit'])
    rep2_300ng_FI6v3_fit    = float(fit_dict[mut]['rep2_300ng_FI6v3_fit'])
    rep1_300ng_FI6v3_E      = rep1_300ng_FI6v3_fit/rep1_no_antibody_fit
    rep2_300ng_FI6v3_E      = rep2_300ng_FI6v3_fit/rep2_no_antibody_fit
    fit_no_antibody = (rep1_no_antibody_fit+rep2_no_antibody_fit)/2 if Input_count >= inputcutoff else -1
    fit_300ng_FI6v3 = (rep1_300ng_FI6v3_fit+rep2_300ng_FI6v3_fit)/2 if Input_count >= inputcutoff else -1
    resist_300ng_FI6v3 = fit_300ng_FI6v3/fit_no_antibody if fit_no_antibody >= fitcutoff and Input_count >= inputcutoff else -1
    if Input_count < inputcutoff: continue
    if mut == 'WT': 
      WT_ipt_count = Input_count
    else:
      all_S_muts.remove(mut)
      outfile.write("\t".join(map(str,[mut,resi,aa,Input_count,fit_no_antibody,fit_300ng_FI6v3,resist_300ng_FI6v3]))+"\n")
  for mut in all_S_muts:
    WTaa  = mut[0]
    mutaa = mut[-1]
    if WTaa != mutaa:
      outfile.write("\t".join(map(str,[mut,mut[1:-1],mut[-1],0,-1,-1,-1]))+"\n")
    else: 
      outfile.write("\t".join(map(str,['WT',mut[1:-1],mut[-1],WT_ipt_count,1,1,1]))+"\n")
  outfile.close()

def main():
  SI06_outfile = 'result/Resist_SI06.tsv'
  Mich15_outfile = 'result/Resist_Mich15.tsv'
  SI06_fit_file = 'result/Fitness_SI06.tsv'
  Mich15_fit_file = 'result/Fitness_Mich15.tsv'
  inputcutoff = 20
  fitcutoff   = 0.5
  SI06_WTresi = ['Q42','I45','N46','G47','I48','T49','V52','H111']
  Mich15_WTresi = ['Q42','I45','D46','K47','I48','T49','V52','H111']
  aas    = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_']
  SI06_fit_dict = TsvWithHeader2Hash(SI06_fit_file)
  Mich15_fit_dict = TsvWithHeader2Hash(Mich15_fit_file)
  fit2escape(SI06_fit_dict, SI06_outfile, inputcutoff, SI06_WTresi, aas, fitcutoff)
  fit2escape(Mich15_fit_dict, Mich15_outfile, inputcutoff, Mich15_WTresi, aas, fitcutoff)

if __name__ == "__main__":
  main()
