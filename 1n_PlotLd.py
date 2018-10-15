__author__ = 'jgwall'

import argparse
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import re

debug = False
# matplotlib.rcParams.update({"font.size":'10'})

def main():
    args = parse_args()
    data = pd.read_table(args.infile, nrows=1000) if debug else pd.read_table(args.infile, dtype=str)
    ld, target = parse_ld(data)
    ld.to_csv(args.outprefix + ".txt", index=False)

    # Set up figure
    fig = plt.figure(figsize=(10, 2))
    grid = gridspec.GridSpec(nrows=1, ncols=1, hspace=0.5, wspace=0.2)
    chromlengths = load_chromlengths(args.chromlengths)

    # Plot Manhatten plots
    ax = fig.add_subplot(grid[:,:], ylabel="Linkage Disequilibrium $(R^2$)")
    plot_ld(ax, ld, chromlengths, sparsify=args.sparsify, zoom=args.zoom, target=target)

    fig.savefig(args.outprefix + ".png", dpi=150)
    fig.savefig(args.outprefix + ".svg", dpi=300)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="TASSEL GWAS output file (site file)")
    parser.add_argument("-o", "--outprefix", help="Output graphic file prefix")
    parser.add_argument("-c", "--chromlengths", help="File of chromosome lengths")
    parser.add_argument("-s", "--sparsify", type=float, default=None, help="What fraction of background GWAS hits to plot")
    parser.add_argument("-z", "--zoom", type=int, help="Chromosome to zoom for output plot")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_chromlengths(infile):
    lengths=dict()
    for line in open(infile, "r"):
        chrom, length, cum = line.strip().split('\t')[:3]   # ':2' in case there are additional columns for some reas
        if chrom.isnumeric(): chrom = int(chrom)
        lengths[chrom] = {"length":int(length), "cum":int(cum)}
    return lengths

def parse_ld(data):
    # Determine which location is the SNP in question
    snp1 = [l1 + "_" + p1 for l1, p1 in zip(data['Locus1'], data['Position1'])]
    snp2 = [l2 + "_" + p2 for l2, p2 in zip(data['Locus2'], data['Position2'])]
    sites, counts = np.unique(snp1 + snp2, return_counts=True)
    target = sites[counts>1]
    
    if(len(target)>1): print("WARNING! More than one target site identified!",target)
    target=target[0]
    
    # Get the position needed for the plot
    chrom = [l1 if s1 != target else l2 for l1, l2, s1 in zip(data['Locus1'], data['Locus2'], snp1)]
    pos = [p1 if s1 != target else p2 for p1, p2, s1 in zip(data['Position1'], data['Position2'], snp1)]
    
    # Make DataFrame
    ld=pd.DataFrame({'Chr':chrom})
    ld['Pos'] = pos
    ld['r2'] = data['R^2']
    
    # Format numerically so works with chromlengths
    ld['Chr'] = [None if c.startswith("SUPER") else int(c) for c in ld['Chr']]     # Hack to deal with non-numeric scaffolds
    ld['Pos'] = [int(p) for p in ld['Pos']]
    ld['r2'] = [float(r) for r in ld['r2']]
    target = [int(x) for x in target.split('_')]
    
    # Remove scaffold rows
    ld = ld.loc[~pd.isnull(ld['Chr']), :].copy()
    
    return(ld, target)


def plot_ld(ax, results, chromlengths, sparsify=None, candidates=None, zoom=None, target=None):

    # Reduce least significant hits, always keeping the top 100 and ones with empirical p-value flags
    if sparsify:
        newsize = math.floor(sparsify * len(results))
        np.random.seed(1)
        tokeep = np.random.choice(range(len(results)), size=newsize, replace=False)
        top_100 = np.argsort(np.array(results['r2']))[:100]

        tokeep = sorted(set(tokeep) | set(top_100))
        results = results.iloc[tokeep,:].copy()

    # Remove NAN values
    results = results.loc[~np.isnan(results['r2']),].copy()

    # x-values
    xvals = list()
    for mychrom, mypos in zip(results['Chr'], results['Pos']):
        xvals.append(mypos + chromlengths[mychrom]['cum'])
    results['x'] = xvals

    # Plot LD
    ax.scatter(x=results['x'], y=results['r2'], s=5, color='darkblue', linewidths=0)
    
    # Plot target point
    target_x = chromlengths[target[0]]['cum'] + target[1]
    ax.axvline(target_x, color='darkred', linewidth=0.5, linestyle='dashed')

    # Add chromosome divisions
    for c in chromlengths:
        ax.axvline(chromlengths[c]['cum'], color="lightgray", zorder=-1)
    ax.axvline(chromlengths[c]['cum'] + chromlengths[c]['length'] , color="lightgray", zorder=-1)

    # Adjust x and y limits
    min_x, max_x = min(results['x']), max(results['x'])
    pad = (max_x - min_x)/20
    ax.set_xlim(left=min_x-pad, right=max_x+pad)
    ax.set_ylim(bottom=0)

    # Set title
    title="Linkage with " + str(target[0]) + ":" + str(target[1])
    ax.set_title(title)
    
    # Alter axis lines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Alter tick marks
    ax.tick_params(axis='y', which='both', right='off', labelsize='x-small')
    ax.tick_params(axis='x', which='both', top='off', bottom='off')

    # Plot chromosome names
    cnames, cpos, chromdiv = list(), list(), set()
    for mychrom in sorted(chromlengths.keys()):
        if mychrom in {'UNKNOWN', 'Pt', 'Mt', '0'}: continue  # skip inconsequential ones
        cnames.append("chr" + str(mychrom))
        cpos.append(chromlengths[mychrom]['cum'] + chromlengths[mychrom]['length'] / 2)
    ax.set_xticks(cpos)
    ax.set_xticklabels(cnames, fontsize="small")

   
    # Zoom if requested
    if zoom:
        min_x, max_x = chromlengths[zoom]['cum'], chromlengths[zoom]['cum'] + chromlengths[zoom]['length']
        pad = 0
        ax.set_xlim(left=min_x-pad, right=max_x+pad)



def set_color(p, cutoff, chrom):
    if p < cutoff and chrom %2 == 0: return "blue"
    if p < cutoff and chrom %2 != 0: return "red"
    if chrom % 2 == 0: return "gray"
    if chrom % 2 != 0: return "silver"
    return "purple" # Should never get here

if __name__ == '__main__': main()
