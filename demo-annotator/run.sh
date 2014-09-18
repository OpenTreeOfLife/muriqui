#!/bin/bash
set -x
#python demo-muriqui.py examples/armadillo-annot.json examples/armadillo-canids.tre
python demo-muriqui.py examples/armadillo-annot.json --tree-file=examples/canids.tre --out-table=examples/canids-out-table.tsv --out-tree=examples/canids-out-tree.tre

