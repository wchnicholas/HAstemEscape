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
  libtype = outfile.rsplit('.')[0].rsplit('_')[1]
  WT_ipt_count = ''
  outfile = open(outfile,'w')
  if libtype == 'S': outfile.write("\t".join(['mut','resi','aa','Input_count','fit_no_antibody',
                                              'fit_2ug_CR9114','fit_10ug_CR9114',
                                              'fit_300ng_FI6v3','fit_2500ng_FI6v3',
                                              'resist_2ug_CR9114','resist_10ug_CR9114',
                                              'resist_300ng_FI6v3','resist_2500ng_FI6v3'])+"\n")
  if libtype == 'D': outfile.write("\t".join(['mut','mut1','mut2','Input_count','fit_no_antibody',
                                              'fit_2ug_CR9114','fit_10ug_CR9114',
                                              'fit_300ng_FI6v3','fit_2500ng_FI6v3',
                                              'resist_2ug_CR9114','resist_10ug_CR9114',
                                              'resist_300ng_FI6v3','resist_2500ng_FI6v3'])+"\n")
  for mut in fit_dict.keys():
    resi  = mut[1:-1]
    aa    = mut[-1] 
    Input_count = float(fit_dict[mut]['Input_count'])
    rep1_no_antibody_fit    = float(fit_dict[mut]['rep1_no_antibody_fit'])
    rep2_no_antibody_fit    = float(fit_dict[mut]['rep2_no_antibody_fit'])
    rep1_2ug_CR9114_fit     = float(fit_dict[mut]['rep1_2ug_CR9114_fit'])
    rep2_2ug_CR9114_fit     = float(fit_dict[mut]['rep2_2ug_CR9114_fit'])
    rep1_10ug_CR9114_fit    = float(fit_dict[mut]['rep1_10ug_CR9114_fit'])
    rep2_10ug_CR9114_fit    = float(fit_dict[mut]['rep2_10ug_CR9114_fit'])
    rep1_300ng_FI6v3_fit    = float(fit_dict[mut]['rep1_300ng_FI6v3_fit'])
    rep2_300ng_FI6v3_fit    = float(fit_dict[mut]['rep2_300ng_FI6v3_fit'])
    rep1_2500ng_FI6v3_fit   = float(fit_dict[mut]['rep1_2500ng_FI6v3_fit'])
    rep2_2500ng_FI6v3_fit   = float(fit_dict[mut]['rep2_2500ng_FI6v3_fit'])
    rep1_2ug_CR9114_E       = rep1_2ug_CR9114_fit/rep1_no_antibody_fit
    rep1_10ug_CR9114_E      = rep1_10ug_CR9114_fit/rep1_no_antibody_fit
    rep2_2ug_CR9114_E       = rep2_2ug_CR9114_fit/rep2_no_antibody_fit
    rep2_10ug_CR9114_E      = rep2_10ug_CR9114_fit/rep2_no_antibody_fit
    rep1_300ng_FI6v3_E      = rep1_300ng_FI6v3_fit/rep1_no_antibody_fit
    rep1_2500ng_FI6v3_E     = rep1_2500ng_FI6v3_fit/rep1_no_antibody_fit
    rep2_300ng_FI6v3_E      = rep2_300ng_FI6v3_fit/rep2_no_antibody_fit
    rep2_2500ng_FI6v3_E     = rep2_2500ng_FI6v3_fit/rep2_no_antibody_fit
    fit_no_antibody = (rep1_no_antibody_fit+rep2_no_antibody_fit)/2 if Input_count >= inputcutoff else -1
    fit_2ug_CR9114 = (rep1_2ug_CR9114_fit+rep2_2ug_CR9114_fit)/2 if Input_count >= inputcutoff else -1
    fit_10ug_CR9114 = (rep1_10ug_CR9114_fit+rep2_10ug_CR9114_fit)/2 if Input_count >= inputcutoff else -1
    fit_300ng_FI6v3 = (rep1_300ng_FI6v3_fit+rep2_300ng_FI6v3_fit)/2 if Input_count >= inputcutoff else -1
    fit_2500ng_FI6v3 = (rep1_2500ng_FI6v3_fit+rep2_2500ng_FI6v3_fit)/2 if Input_count >= inputcutoff else -1
    resist_2ug_CR9114 = fit_2ug_CR9114/fit_no_antibody if fit_no_antibody >= fitcutoff and Input_count >= inputcutoff else -1
    resist_10ug_CR9114 = fit_10ug_CR9114/fit_no_antibody if fit_no_antibody >= fitcutoff and Input_count >= inputcutoff else -1
    resist_300ng_FI6v3 = fit_300ng_FI6v3/fit_no_antibody if fit_no_antibody >= fitcutoff and Input_count >= inputcutoff else -1
    resist_2500ng_FI6v3 = fit_2500ng_FI6v3/fit_no_antibody if fit_no_antibody >= fitcutoff and Input_count >= inputcutoff else -1
    if Input_count < inputcutoff: continue
    if mut == 'WT': 
      WT_ipt_count = Input_count
    else:
      if libtype == 'S':
        all_S_muts.remove(mut)
        outfile.write("\t".join(map(str,[mut,resi,aa,Input_count,fit_no_antibody,fit_2ug_CR9114,fit_10ug_CR9114,fit_300ng_FI6v3,fit_2500ng_FI6v3,resist_2ug_CR9114,resist_10ug_CR9114,resist_300ng_FI6v3,resist_2500ng_FI6v3]))+"\n")
      elif libtype == 'D':
        if '/' not in mut: continue
        mut1 = mut.rsplit('/')[0]
        mut2 = mut.rsplit('/')[1]
        outfile.write("\t".join(map(str,[mut,mut1,mut2,Input_count,fit_no_antibody,fit_2ug_CR9114,fit_10ug_CR9114,fit_300ng_FI6v3,fit_2500ng_FI6v3,resist_2ug_CR9114,resist_10ug_CR9114,resist_300ng_FI6v3,resist_2500ng_FI6v3]))+"\n")
      else: print "Something is wrong with libtype"; sys.exit()
  if libtype == 'S': 
    for mut in all_S_muts:
      WTaa  = mut[0]
      mutaa = mut[-1]
      if WTaa != mutaa:
        outfile.write("\t".join(map(str,[mut,mut[1:-1],mut[-1],0,-1,-1,-1,-1,-1,-1,-1,-1,-1]))+"\n")
      else: 
        outfile.write("\t".join(map(str,['WT',mut[1:-1],mut[-1],WT_ipt_count,1,1,1,1,1,1,1,1,1]))+"\n")
  outfile.close()

def main():
  S_outfile = 'result/Resist_S.tsv'
  D_outfile = 'result/Resist_D.tsv'
  Sfit_file = 'result/Fitness_S.tsv'
  Dfit_file = 'result/Fitness_D.tsv'
  inputcutoff = 20
  fitcutoff   = 0.5
  WTresi = ['Q42','I45','D46','Q47','I48','N49','L52','T111']
  aas    = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_']
  Sfit_dict = TsvWithHeader2Hash(Sfit_file)
  Dfit_dict = TsvWithHeader2Hash(Dfit_file)
  fit2escape(Sfit_dict, S_outfile, inputcutoff, WTresi, aas, fitcutoff)
  fit2escape(Dfit_dict, D_outfile, inputcutoff, WTresi, aas, fitcutoff)

if __name__ == "__main__":
  main()
