# Write matplotlib colormaps to ascii file.
# Retrieve colormap names with argument 'listCMapNames'
# Version 1.1; K. Schumacher 08.2018

import matplotlib.pyplot as plt
import numpy as np
import argparse


def get_cmap(name, nVal):
    cmp = plt.get_cmap(name, nVal)
    return cmp(range(nVal))

def write_cmap(file, cmap):
    np.savetxt(file, cmap)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmName', type=str, choices=['listCMapNames']+plt.colormaps())
    ap.add_argument('-o', type=str, dest='outfile')
    ap.add_argument('-n', type=int, default=256, dest='nVal')
    
    args = ap.parse_args()
    
    if args.outfile is None:
        args.outfile = args.cmName + '.txt'

    if args.cmName in 'listCMapNames':
        with open(args.outfile,'w') as f:
            f.writelines(c+'\n' for c in plt.colormaps())
    else:
        cmp = get_cmap(args.cmName, args.nVal)
        write_cmap(args.outfile, cmp)
    