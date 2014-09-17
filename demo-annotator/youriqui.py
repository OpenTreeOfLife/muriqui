#!/usr/bin/env python
import dendropy
import codecs
import json
import sys
try:
    tree_file, annotations_file = sys.argv[1:]
except:
    sys.exit('expecting 2 arguments: a tree file and a JSON file of annotations')
tree_list = dendropy.TreeList.get_from_path(tree_file,
                                            'newick')
for tree in tree_list:
    tree.print_plot()
    tree.encode_splits(delete_outdegree_one=False)
    for node in tree.preorder_node_iter():
        print node.edge.split_bitmask

a_f = codecs.open(annotations_file, 'rU', encoding='utf-8')
annot_list = json.load(a_f)
if not isinstance(annot_list, list):
    annot_list = [annot_list]
for annotation in annot_list:
    print annotation
