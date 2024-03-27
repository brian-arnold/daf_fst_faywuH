#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import numpy as np
import seaborn as sns

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

  pop_set = "CEZ_BLR_YLS"
  dafs = []
  positions = []
  ddaf_file = f"../{pop_set}_WIN21_OFF1_MINAF1_DDAF.txt"
  with open(ddaf_file, 'r') as d:
    for line in d:
      if line.startswith("chrX"):
        
        line = line.strip().split()
        pos = int(line[1])
        daf = float(line[5])
        dafs.append(daf)
        positions.append(pos)

  winsize = 50
  offset = 10

  windows = sliding_windows(dafs, winsize, offset)
  dafs_means = []
  for win in windows:
    dafs_means.append(np.mean(win))

  windows = sliding_windows(positions, winsize, offset)
  positions_means = []
  for win in windows:
    positions_means.append(np.median(win))


  print(len(dafs_means))
  print(len(positions_means))

  # make a scatter plot in seaborn where the x values come from positions_means and the y values come from dafs_means, and save it to the current directory
  scat = sns.scatterplot(x=positions_means, y=dafs_means, linewidth=0, alpha=0.5)
  sns.despine()
  # save figure
  fig = scat.get_figure()
  fig.savefig(f"daf_vs_position_{pop_set}.png", dpi=300)


  



if __name__ == '__main__':
  main()
