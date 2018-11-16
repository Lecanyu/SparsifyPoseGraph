#!/usr/bin/env python

import sys, os.path

datasets = {
    "duderstadt": ("se3", 544),
    "eecs": ("se3", 614),
    "flat_duderstadt": ("se2", 544),
    "flat_eecs": ("se2", 614),
    "intel": ("se2", 942),
    "killian": ("se2", 5488),
    "manhattan": ("se2", 3499),
    "parking_fixed": ("se3", 1660),
    "sphere": ("se3", 2499)
}

def line2filename(line, resultsdir):
    tokens = line.split(" ")
    dataset = tokens[1].split(".")[0].split("/")[1]
    type = tokens[0] if tokens[0] == "glc" else datasets[dataset][0]
    profile = tokens[2]
    topology = tokens[3]
    linpoint = "l" if tokens[4] == "local" else "g"
    sparsity = tokens[5]
    
    return "%s/%s/%s/%s/%s_%s_%s.kld" % ( \
        resultsdir, profile, sparsity, dataset, type, topology, linpoint)

def lastline(fname):
    dataset = fname.split("/")[3]
    lines = open(fname).readlines()
    if len(lines) > 0:
        line = lines[-1].strip()
    else:
        line = ""
    if len(line) > 0:
        tokens = line.split(" ")
        return (dataset, int(tokens[0]), float(tokens[1]))
    else:
        return (dataset, 0, 0)

def main(argv):
    if len(argv) >= 2:
        if len(argv) == 3:
            resultsdir = argv[2]
        else:
            resultsdir = "full_results"
        for line in open(argv[1]):
            line = line.strip()
            if line.startswith("#") or len(line) == 0:
                continue
            
            fname = line2filename(line, resultsdir)
            if os.path.isfile(fname):
                dataset, iter, value = lastline(fname)
                if iter < datasets[dataset][1]:
                    print "%4d/%4d %s" % (iter, datasets[dataset][1], fname)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
