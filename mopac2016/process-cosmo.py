#!/usr/bin/env python3

import sys
import os
import numpy as np

def weightBinSigmas(sigmavals, sigmas_grid):
    """
    """
    # Regular grid, so every bin is of same width
    bin_width = sigmas_grid[1] - sigmas_grid[0]

    psigmaAs = np.zeros_like(sigmas_grid)
    for sigma, area in sigmavals:
        # Check sigma
        if sigma < np.min(sigmas_grid):
          continue
          #  raise ValueError('Sigma [{0:g}] is less than minimum of grid [{1}]'.format(sigma, np.min(sigmas_grid)))
        if sigma > np.max(sigmas_grid):
          continue
          #  raise ValueError('Sigma [{0:g}] is greater than maximum of grid [{1}]'.format(sigma, np.max(sigmas_grid)))
        # The index to the left of the point in sigma
        left = int((sigma-sigmas_grid[0])/bin_width)
        # Weighted distance from the left edge of the cell to the right side
        w_left = (sigmas_grid[left+1]-sigma)/bin_width
        # Add the area into the left and right nodes of the bin, each part weighted by the value
        # of sigma, if equal to the left edge of the bin, then the w_left=1, if the right side,
        # w_left = 0
        psigmaAs[left] += area*w_left
        psigmaAs[left+1] += area*(1.0-w_left)
    return psigmaAs

# Main

# Read the command line input
fileName = "FOR005.cos"
smiles = "input"
bin_width = float(sys.argv[1])

# Read and parse the input file to get list of sigma values

try:
  inputFile = open(fileName, "r").read()
  inputLines = inputFile.splitlines()

  # Extract the sigma data from mopac .cos output
  sigmas = []
  totalVolume = 0
  totalArea = 0
  for count, line in enumerate(inputLines):
    if line.startswith('          COSMO AREA'):
      fields = list(filter(None, line.split(' ')))
      totalArea = fields[3]
    if line.startswith('          COSMO VOLUME'):
      fields = list(filter(None, line.split(' ')))
      totalVolume = fields[3]
    if line.startswith('           SEGMENT DATA:'):
      startIndex = count + 2
      for index in range(startIndex, len(inputLines)):
        values = list(filter(None, inputLines[index].split(' ')))
        floatValues = [float(value) for value in values]
        sigma = floatValues[8]
        area = floatValues[7]
        pair = [floatValues[8], floatValues[7]]
        sigmas.append(pair)
except:
  print('No Such File')


sigmaRanges = np.arange(-0.025, 0.025+0.0001, bin_width) # [e/A^2]

# Weight and Bin the Sigmas over Specified bins
psigmaAs = weightBinSigmas(sigmas, sigmaRanges)

# Create the output sigma profile
out = '# File Name: %s\n# SMILES: %s\n# Total Volume: %s\n# Total Area: %s\n' % (fileName, smiles, totalVolume, totalArea)
for index in range(len(psigmaAs)):
  out += '{0:0.6f} {1:0.6f}\n'.format(sigmaRanges[index], psigmaAs[index])

outputFileName = "sigma.txt"

with open(outputFileName, 'w') as outputFile:
  outputFile.write(out)

