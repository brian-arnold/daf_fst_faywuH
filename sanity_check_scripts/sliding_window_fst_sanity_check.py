#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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

def main():

  positions = []
  fst = defaultdict(list)
  fst_file = "../CEZ_BLR_YLS_WIN21_OFF1_MINAF1_FST.txt"
  with open(fst_file, 'r') as f:
    for line in f:
      if line.startswith("chrX"):

        line = line.strip().split()
        midpos = int(line[4])
        try:
          fst['BLR_CEZ'].append(float(line[6]))
        except ValueError:
          fst['BLR_CEZ'].append(np.nan)

        try:
          fst['BLR_YLS'].append(float(line[7]))
        except ValueError:
          fst['BLR_YLS'].append(np.nan)

        try:
          fst['CEZ_YLS'].append(float(line[8]))
        except ValueError:
            fst['CEZ_YLS'].append(np.nan)

        positions.append(midpos)


  # make a scatter plot in seaborn where the x values come from positions_means and the y values come from dafs_means, and save it to the current directory
  fig, axs = plt.subplots(3,1, figsize=(10, 10))
  scat = sns.scatterplot(x=positions, y=fst['BLR_CEZ'], linewidth=0, alpha=0.1, ax=axs[0])
  axs[0].set_title("BLR_CEZ")
  sns.scatterplot(x=positions, y=fst['BLR_YLS'], linewidth=0, alpha=0.1, ax=axs[1])
  axs[1].set_title("BLR_YLS")
  sns.scatterplot(x=positions, y=fst['CEZ_YLS'], linewidth=0, alpha=0.1, ax=axs[2])
  axs[2].set_title("CEZ_YLS")
  sns.despine()
  # save figure
  scat.get_figure().savefig("fst_vs_position.png", dpi=300)
  


  



if __name__ == '__main__':
  main()
