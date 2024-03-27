#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import math
from itertools import combinations
import numpy as np

def compute_segsites(win):
  return len(win)

def compute_thetaW(S, a1):
  return S/a1

def compute_thetaPi(win, n):
  tPi = 0
  for w in win:
      # 2rd element of tuple is AF
      tPi += (w[1]*(n-w[1]))/(n*(n-1)/2)
  return tPi

def compute_thetaL(win, n):
  tL = 0
  for w in win:
    # 2rd element of tuple is AF
    tL += w[1]/(n-1)
  return tL

def compute_faywuh(win, constants):
  n = constants['two_n']
  a1, a2, b2 = constants['a1'], constants['a2'], constants['b2']
  S = compute_segsites(win)
  tW = compute_thetaW(S, a1)
  tPi = compute_thetaPi(win, n)
  tL = compute_thetaL(win, n)

  thetaSquared = S*(S-1)/((a1**2) + a2)
  var_pt1 = (n-2)*tW/(6*(n-1))
  var_pt2 = (18*(n**2))*((3*n) + 2)*b2 - ((88*(n**3))+(9*(n**2))-(13*n)+6)
  var_pt3 = thetaSquared/(9*n*((n-1)**2))
  var = var_pt1 + (var_pt2*var_pt3)
  if math.sqrt(var) > 0:
    H = (tPi - tL)/math.sqrt(var)
    return(H)
  else:
    return("NA")
  
def polarize_allele_frequency(index2Group, variant, two_n, two_n_OG):

  AF = defaultdict(int)
  for i in index2Group:	
    grp = index2Group[i]
    AF[grp] += variant.gt_types[i]

  ancAllele = ""
  # if reference population has more of an ALT allele, use this as ancestral
  if AF['outgroup'] > (two_n_OG/2):
    ancAllele = str(variant.ALT[0])
  else:
    ancAllele = variant.REF

  if ancAllele == variant.REF:
    return AF['ingroup']
  else:
    return two_n - AF['ingroup']

def make_constants(n):
  a1 = 0
  a2 = 0
  b2 = 0
  e1 = 0
  e2 = 0
  # range goes until end-1, so use n for end bc you want n-1
  for i in range(1, n):
    a1 += 1/i
    a2 += 1/i**2
  # used in FayWuH variance calculation, Zeng et al 2006 Genetics
  # b2 in bn+1, where n is +1 other constants
  for i in range(1, n+1):
    b2 += 1/i**2

  # For Tajima's D
  e1 = (1/a1)*((n+1)/(3*(n-1)) - (1/a1))
  e2 = (1/(a1**2 + a2))*((2*((n**2) + n + 3))/(9*n*(n-1)) - (n+2)/(n*a1) + (a2)/(a1**2))

  constants = defaultdict()
  constants['a1'] = a1
  constants['a2'] = a2
  constants['b2'] = b2  
  constants['e1'] = e1
  constants['e2'] = e2
  constants['two_n'] = n
  return constants


def compute_hudson_fst(win, two_n):
  # win is a list of variants
  fst_pop_pair = defaultdict(float)
  snp_pop_pair = defaultdict(float)

  pop_pairs = list(combinations(sorted(two_n.keys()), 2))
  for pair in pop_pairs:
    pop1, pop2 = pair
    numerators, denominators = [], []
    for w in win:
      # 2nd element of w is the dict of allele counts for all populations
      p1 = w[1][pop1]/two_n[pop1]
      p2 = w[1][pop2]/two_n[pop2]

      # if site invariant for this population pair, FST isn't defined
      if p1+p2 == 0 or p1+p2 == 2:
        continue

      # computing Hudson's FST from Bhatia, Patterson, Sankararaman, Price (2013), Genome Research
      # https://genome.cshlp.org/content/23/9/1514.full
      numerator = (p1 - p2)**2 - ((p1*(1-p1))/(two_n[pop1]-1)) - ((p2*(1-p2))/(two_n[pop2]-1))
      denominator = p1*(1-p2) + p2*(1-p1)
      numerators.append(numerator)
      denominators.append(denominator)
    # Bhatia et al 2013 Genome Research recommend using a "ratio of averages" in which numerator and denomiator are averaged separately
    if len(numerators) > 0:
      assert len(numerators) == len(denominators)
      fst_pop_pair["_".join(pair)] = np.sum(numerators)/np.sum(denominators)
      snp_pop_pair["_".join(pair)] = len(numerators)
    else:
      fst_pop_pair["_".join(pair)] = "NA"
      snp_pop_pair["_".join(pair)] = 0 
  return fst_pop_pair, snp_pop_pair

def compute_DDAF(AFs, two_n, args):
  focal_pop_AF = AFs[args.focal_pop]/two_n[args.focal_pop]
  sister_pop_AF = AFs[args.sister_pop]/two_n[args.sister_pop]
  outgroup_pop_AF = AFs[args.outgroup_pop]/two_n[args.outgroup_pop]

  return focal_pop_AF - ((sister_pop_AF+outgroup_pop_AF)/2)
