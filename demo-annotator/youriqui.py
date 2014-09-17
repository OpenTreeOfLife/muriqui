#!/usr/bin/env python
import dendropy
import codecs
import json
import sys
import os
SCRIPT_NAME = os.path.split(sys.argv[0])[1]
class Reason(object):
    NO_INC_DESIGNATORS_IN_TREE = 0
    SUCCESS = 1
    MRCA_HAS_EXCLUDED = 2

def debug(msg):
    sys.stderr.write('{s}: {m}\n'.format(s=SCRIPT_NAME, m=msg))
def taxa_in_tree(tree, taxa_list):
    t = []
    for i in taxa_list:
        ind = tree.label2index.get(i)
        if ind is not None:
            t.append(tree.taxon_namespace[ind])
    return t
def bits_in_tree(tree, taxa_list):
    t = taxa_in_tree(tree, taxa_list)
    return [tree.label2bit[i.label] for i in t]
def add_node_based_phyloreferenced_annotation(tree, annotation):
    return False, -1
def add_stem_based_phylorefenced_annotation(tree, annotation):
    in_tree = taxa_in_tree(tree, annotation.des)
    exclude = bits_in_tree(tree, annotation.exclude_ancs_of)
    if not exclude:
        return tree.seed_node, Reason.SUCCESS
    exc_bit_set = 0
    for e in exclude:
        exc_bit_set |= e
    if len(in_tree) == 0:
        return False, Reason.NO_INC_DESIGNATORS_IN_TREE
    if len(in_tree) == 1:
        mrca = tree.find_node_with_taxon_label(in_tree[0].label)
        assert mrca is not None
    else:
        mrca = tree.mrca(taxa=in_tree)
    if mrca.edge.split_bitmask & exc_bit_set:
        return False, Reason.MRCA_HAS_EXCLUDED
    deepest_valid = mrca
    curr = mrca.parent_node
    while curr is not None and not (curr.edge.split_bitmask & exc_bit_set):
        deepest_valid = curr
        curr = curr.parent_node
    return curr.edge, Reason.SUCCESS
def add_phyloreferenced_annotation(tree, annotation):
    if annotation.exclude_ancs_of:
        return add_stem_based_phylorefenced_annotation(tree, annotation)
    return add_node_based_phylorefenced_annotation(tree, annotation)
class PhyloReferencedAnnotation(object):
    def __init__(self, serialized):
        self.des = [str(i) for i in serialized['descendants']]
        self.exclude_ancs_of = [str(i) for i in serialized.get('excludes_ancestors_of', [])]
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
        tree.label2index = {}
        tree.label2bit = {}
        curr_bit = 1
        for n, taxon in enumerate(tree.taxon_namespace):
            assert taxon.label not in tree.label2index
            tree.label2index[taxon.label] = n
            tree.label2bit[taxon.label] = curr_bit
            curr_bit <<= 1

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
            x, reason = add_phyloreferenced_annotation(tree, a)
            if x:
                num_added += 1
            else:
                debug('Annotation {a} could not be added to tree {t}'.format(a=annot_index, t=tree_index))
            num_tried += 1
        debug('{a}/{t} annotations added to tree {i}'.format(a=num_added, t=num_tried, i=tree_index))
