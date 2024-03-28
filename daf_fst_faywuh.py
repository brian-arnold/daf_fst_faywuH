#!/usr/bin/python -tt

import re
import sys
import os
from cyvcf2 import VCF
import math
from collections import defaultdict
import numpy as np
#sys.path.append(os.path.join(os.path.dirname(__file__)))
#sys.path.append('./')
from parse_args import args
from itertools import combinations
import popgen_stats

def filterVariantDepthMissData(gt_depths, gt_types, index2Group, Groups, Filter):
	PASS = 1
	# USE LIST; YOU NEED INDICES OF MISSING INDIVIDUALS LATER SO YOU DON'T USE THEM FOR CALCULATING AF 
	sampleSizes_missing = defaultdict(list)

	for i in index2Group:
		if (gt_depths[i] < Filter['MinDepthPerInd']) or (gt_types[i] == 3): # has low depth or missing genotype (unlikely for reasonable depth thresholds), but for which category?
			grp = index2Group[i]
			sampleSizes_missing[grp].append(i)
	for grp in sampleSizes_missing:
		if len(sampleSizes_missing[grp]) > Filter['MaxMissIndPerCategory']: 
			PASS = 0
	return PASS, sampleSizes_missing

def sliding_windows(iterable, size, offset):
  it = iter(iterable)
  win = []
  # create first window with 'size' elements
  for i in range(0, size):
    try:
      win.append(next(it))
    except StopIteration:
      return
  yield win
  # create rest of windows
  # at each iteration, remove first element of window and add next element to window
  # once j reaches offset, yield window
  j = 0
  for i in it:
    win = win[1:] + [i]
    j += 1
    if j % offset == 0:
      yield win
   
def get_positions(win):
  positions = np.array([w[0] for w in win])
  return np.min(positions), np.max(positions), np.mean(positions)

def get_AFs(index2Group, variant, two_n, args):
  # returns a dict of allele frequencies for each population
  AF = defaultdict(int)

  for i in index2Group:	
    grp = index2Group[i]
    AF[grp] += variant.gt_types[i]

  # polarize allele frequencies using the major allele in the outgroup population
  # if ALT frequency greater than 0.5, invert AFs
  if AF[args.outgroup_pop] > (two_n[args.outgroup_pop]/2):
    for pop in AF:
      new = two_n[pop] - AF[pop]
      AF[pop] = new

  return AF



def parse_pop_file(pop_file):
  groups = defaultdict(list)
  for line in open(pop_file, "r"):
    line = line.strip()
    if line == "":
      continue
    line = line.split()
    groups[line[1]].append(line[0])
  return groups

def main():

  assert args.min_AF > 0, "min_AF parameter must be greater than 0, otherwise invariant sites get included"
  assert  args.winsize % 2 == 1, "winsize must be odd"
  vcf = VCF(args.vcf, gts012=True) # gts012=True makes value of 3 UNKNOWN, with 0,1,2 corresponding to numb ALT alleles

  groups = parse_pop_file(args.pop_file)
  print(groups)
  assert args.focal_pop in groups, "focal population not in pop file!"
  assert args.sister_pop in groups, "sister population not in pop file!"

  all_samples = [ind for pop in groups for ind in groups[pop] ]
  print(all_samples)

  vcf.set_samples( all_samples )

  n = {pop : len(groups[pop]) for pop in groups}
  two_n = {pop : len(groups[pop])*2 for pop in groups}
  pop_pair_two_n = two_n[args.focal_pop] + two_n[args.sister_pop]
  print(n)
  print(two_n, sum([i for i in two_n.values()]))
  
  #############
  # for each index in vcf.samples list, which group does it correspond to
  #############
  index2Group = {} 
  for index, samp in enumerate(vcf.samples):
    for grp in groups:
      if samp in groups[grp]:
        index2Group[index] = grp
  print('index2Group: ', index2Group)
  # make Dict to get idea of missing individuals per group, sums all individuals w/ missing data, across all sites
  filteredSamplesByGrp = defaultdict(int)
  FilterDP = defaultdict(int)
  FilterDP['MinDepthPerInd'] = 4

  # polarized_AFs[chromosome] = list of 2-tuples, position and polarized AF
  focal_pop_AFs = defaultdict(list) # keys are chromosomes, values are lists of tuples (position, AF for focal pop)
  sisters_pop_AFs = defaultdict(list) # keys are chromosomes, values are lists of tuples (position, AF dict for all pops, effect, gene name)
  for variant in vcf:
    # only look at biAllelic Sites
    if len(variant.ALT) == 1 and variant.var_type == "snp":
      # go thru gt_depths, check if greater than Filter, tally up miss data for each category
      PASS_depthMissData, sampleSizes_missing = filterVariantDepthMissData(variant.gt_depths, variant.gt_types, index2Group, groups, FilterDP)
      # gt_types is array of 0,1,2,3==HOM_REF, HET, HOM_ALT, UNKNOWN 
      if PASS_depthMissData != 1:
        continue

      AFs = get_AFs(index2Group, variant, two_n, args)
      # check if site variable across both pops
      total_AF = sum(AFs.values())
      assert total_AF >=0 and total_AF <= sum([i for i in two_n.values()])
      # collect AFs for positions variable in focal pop
      if AFs[args.focal_pop] >= args.min_AF and AFs[args.focal_pop] <= (two_n[args.focal_pop] - args.min_AF):
        # STORE ALLELE FREQS AS DICT E.G. have tuple with (position, {pop1: AF, pop2: AF})
        # assures populations never get mixed up
        info = tuple( [variant.end, AFs[args.focal_pop]] )
        focal_pop_AFs[variant.CHROM].append( info )

      pop_pair_AFs = AFs[args.focal_pop] + AFs[args.sister_pop]
      if pop_pair_AFs >= args.min_AF and pop_pair_AFs <= (pop_pair_two_n - args.min_AF):
        effect = "NA"
        geneName = "NA"
        for field in variant.INFO:
          if field[0] == 'ANN':
            info = field[1]
            infoList = info.split('|')
            effect = infoList[1]
            geneName = infoList[3]
        info = tuple( [variant.end, AFs, effect, geneName] )
        sisters_pop_AFs[variant.CHROM].append( info )

  # compute FST statistic for each window
  results = []
  for chrom in sisters_pop_AFs:
    windows = sliding_windows(sisters_pop_AFs[chrom], args.winsize, args.offset)
    for win in windows:
      # give FST function total number of chromosomes samples in each population
      # allele frequencies are currently integers, so divide by 2N
      fst_pop_pair_dict, snp_pop_pair_dict = popgen_stats.compute_hudson_fst(win, two_n)
      minPos, maxPos, midPos = get_positions(win)
      midSNPpos = win[ int(len(win)/2) ][0] # don't need to add 1 since index starts at 0
      assert midSNPpos > minPos and midSNPpos < maxPos, "midSNPpos not in window!"
      assert len(win) == snp_pop_pair_dict["_".join(sorted([args.focal_pop, args.sister_pop]))], "number of snps in window not equal to expected!"
      snps = [snp_pop_pair_dict[ppair] for ppair in sorted(snp_pop_pair_dict)]
      fsts = [fst_pop_pair_dict[ppair] for ppair in sorted(fst_pop_pair_dict)]
      # use python map function to convert all variables to strings
      r = [chrom, minPos, maxPos, midPos, midSNPpos, len(win)] + fsts + snps
      results.append( list(map( str, r )) )
  # write FST results to file
  with open(args.output_fst, 'w') as outfile:
    names_fst = [f'Fst_{ppair}' for ppair in sorted(fst_pop_pair_dict)]
    names_snp = [f'SNP_{ppair}' for ppair in sorted(snp_pop_pair_dict)]
    outfile.write("\t".join(["chrom", "minPos", "maxPos", "midPos", "midSNPpos", "winLen"] + names_fst + names_snp + ["\n"]))
    outfile.write("\n".join(["\t".join(r) for r in results]) + "\n")
  
  # compute DAF
  results = []
  for chrom in sisters_pop_AFs:
    for snp in sisters_pop_AFs[chrom]:
      DDAF = popgen_stats.compute_DDAF(snp[1], two_n, args)
      derived_AFs = [snp[1][pop] for pop in sorted(snp[1])]
      effect, geneName = snp[2], snp[3]
      r = [chrom, snp[0]] + derived_AFs + [DDAF, effect, geneName]
      results.append( list(map( str, r )) )
  with open(args.output_ddaf, 'w') as outfile:
    names = sorted(groups)
    outfile.write("\t".join(["chrom", "pos"] + names + [ "ddaf", "\n"]))
    outfile.write("\n".join(["\t".join(r) for r in results]) + "\n")     

  # compute Fay and Wu's H statistic for each window
  constants = popgen_stats.make_constants(two_n[args.focal_pop])
  results = []
  for chrom in focal_pop_AFs:
    windows = sliding_windows(focal_pop_AFs[chrom], args.winsize, args.offset)
    for win in windows:
      faywuh = popgen_stats.compute_faywuh(win, constants)
      minPos, maxPos, midPos = get_positions(win)
      midSNPpos = win[ int(len(win)/2) ][0] # don't need to add 1 since index starts at 0
      # use python map function to convert all variables to strings
      r = [chrom, minPos, maxPos, midPos, midSNPpos, len(win), faywuh]
      results.append( list(map(str, r)) )
  # write Fay and Wu's H results to file
  with open(args.output_faywuh, 'w') as outfile:
    outfile.write("\t".join(["chrom", "minPos", "maxPos", "midPos", "midSNPpos", "winLen", "FayWuH", "\n"]))
    outfile.write("\n".join(["\t".join(r) for r in results]) + "\n")




  """
  # write to file allele frequencies per population, per site
  with open(args.output_AF, 'w') as outfile:
    pops = sorted(focal_pop_AFs.keys())
    outfile.write("\t".join(["chrom", "pos", "AF_pop1", "AF_pop2","\n"]))
    for chrom in pop_AFs:
      for pos, AF1, AF2 in pop_AFs[chrom]:
        outfile.write("\t".join([str(chrom), str(pos), str(AF1/two_n_p1), str(AF2/two_n_p2)]) + "\n")
  """
if __name__ == '__main__':
  main()
