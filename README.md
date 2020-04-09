## Getting Started

```sh
## If you can run paftools.js from minimap2, you already have k8 installed. If not:

## install k8 without conda:
# curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
# cp k8-0.2.4/k8-`uname -s` k8  # or copy it to a directory on your $PATH

## install k8 via conda:
# conda install minimap2        # k8 comes with minimap2

k8 calN50.js ctg.fa             # compute auN and N50 from FASTA
k8 calN50.js ctg.fa.fai         # faidx index (only first two fields are needed)
calN50.js seq.gfa.gz            # if k8 and calN50.js are on $PATH
calN50.js -L3.1g seq.lens       # compute auNG and NG50 for a 3.1Gbp genome
calN50.js -f ref.fa.fai seq.fa  # or get the genome size from a .fai file
```

## Introduction

calN50.js is a simple script to calculate N50/NG50 and [auN][auN]/auNG.

[auN]: http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity
