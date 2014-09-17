#!/usr/bin/env python
import sys
try:
    tree_file, annotations_file = sys.argv[1:]
except:
    sys.exit('expecting 2 arguments: a tree file and a JSON file of annotations')




