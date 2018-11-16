#!/usr/bin/env python

import sys, random

def main(argv):
    if len(argv) >= 3:
        lines = [ s for s in open(argv[1]).readlines() if not s.startswith("#") ]
        random.shuffle(lines)
        n = int(argv[2])
        sizes = [len(lines) / n] * n
        for i in range(len(lines) % n):
            sizes[i] += 1
        k = 0
        for i in range(len(sizes)):
            open("input-%d.txt" % i, "w").writelines(lines[k:k+sizes[i]])
            k += sizes[i]
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
