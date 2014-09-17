#!/usr/bin/env python
import dendropy
import codecs
import json
import sys
import os
SCRIPT_NAME = os.path.split(sys.argv[0])[1]
def debug(msg):
    sys.stderr.write('{s}: {m}\n'.format(s=SCRIPT_NAME, m=msg))
def add_phyloreferenced_annotation(tree, annotation):
    return False
class PhyloReferencedAnnotation(object):
    def __init__(self, serialized):
        self.des = serialized['descendants']
        self.exclude_ancs_of = serialized.get('excludes_ancestors_of', [])
if __name__ == '__main__':
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

    for tree_index, tree in enumerate(tree_list):
        num_tried = 0
        num_added = 0
        for annot_index, annotation in enumerate(annot_list):
            a = PhyloReferencedAnnotation(annotation)
            x = add_phyloreferenced_annotation(tree, a)
            if x:
                num_added += 1
            else:
                debug('Annotation {a} could not be added to tree {t}'.format(a=annot_index, t=tree_index))
            num_tried += 1
        debug('{a}/{t} annotations added to tree {i}'.format(a=num_added, t=num_tried, i=tree_index))
