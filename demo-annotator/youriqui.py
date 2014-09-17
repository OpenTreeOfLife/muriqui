#!/usr/bin/env python
import dendropy
from dendropy.utility import container
import codecs
import json
import sys
import os
SCRIPT_NAME = os.path.split(sys.argv[0])[1]
class Reason(object):
    NO_INC_DESIGNATORS_IN_TREE = 0
    SUCCESS = 1
    MRCA_HAS_EXCLUDED = 2
def mod_encode_splits(tree, create_dict=True, delete_outdegree_one=True, internal_node_taxa=False):
    """
    Processes splits on a tree, encoding them as bitmask on each edge.
    Adds the following to each edge:
        - `split_bitmask` : a rooted split representation, i.e. a long/bitmask
            where bits corresponding to indices of taxa descended from this
            edge are turned on
    If `create_dict` is True, then the following is added to the tree:
        - `split_edges`:
            [if `tree.is_rooted`]: a dictionary where keys are the
            splits and values are edges.
            [otherwise]: a container.NormalizedBitmaskDictionary where the keys are the
            normalized (unrooted) split representations and the values
            are edges. A normalized split_mask is where the split_bitmask
            is complemented if the right-most bit is not '0' (or just
            the split_bitmask otherwise).
    If `delete_outdegree_one` is True then nodes with one
        will be deleted as they are encountered (this is required
        if the split_edges dictionary is to refer to all edges in the tree).
        Note this will mean that an unrooted tree like '(A,(B,C))' will
        be changed to '(A,B,C)' after this operation!
    """
    taxon_namespace = tree.taxon_namespace
    if taxon_namespace is None:
        taxon_namespace = tree.infer_taxa()
    if create_dict:
        tree.split_edges = {}
        split_map = tree.split_edges
        # if tree.is_rooted:
        #     tree.split_edges = {}
        # else:
        #     atb = taxon_namespace.all_taxa_bitmask()
        #     d = container.NormalizedBitmaskDict(mask=atb)
        #     tree.split_edges = d
        # split_map = tree.split_edges
    if not tree.seed_node:
        return

    if delete_outdegree_one:
        sn = tree.seed_node
        if not tree.is_rooted:
            if len(sn._child_nodes) == 2:
                tree.deroot()
        while len(sn._child_nodes) == 1:
            c = sn._child_nodes[0]
            if len(c._child_nodes) == 0:
                break
            try:
                sn.edge.length += c.edge.length
            except:
                pass
            sn.remove_child(c)
            for gc in c._child_nodes:
                sn.add_child(gc)

    for edge in tree.postorder_edge_iter():
        cm = 0
        h = edge.head_node
        child_nodes = h._child_nodes
        nc = len(child_nodes)
        ocm = None
        if nc > 0:
            if nc == 1 and delete_outdegree_one and edge.tail_node:
                p = edge.tail_node
                assert(p)
                c = child_nodes[0]
                try:
                    c.edge.length += edge.length
                except:
                    pass
                pos = p._child_nodes.index(h)
                p.insert_child(pos, c)
                p.remove_child(h)
            else:
                for child in child_nodes:
                    cm |= child.edge.split_bitmask
                if internal_node_taxa:
                    t = edge.head_node.taxon
                    if t:
                        ocm = taxon_namespace.taxon_bitmask(t)
                        cm |= ocm
        else:
            t = edge.head_node.taxon
            if t:
                cm = taxon_namespace.taxon_bitmask(t)
        edge.split_bitmask = cm
        if create_dict:
            split_map[cm] = edge
            if ocm is not None:
                split_map[ocm] = edge
    # create normalized bitmasks, where the full (tree) split mask is *not*
    # all the taxa, but only those found on the tree
    if not tree.is_rooted:
        mask = tree.seed_node.edge.split_bitmask
        d = container.NormalizedBitmaskDict(mask=mask)
        for k, v in tree.split_edges.items():
            d[k] = v
        tree.split_edges = d

def debug(msg):
    sys.stderr.write('{s}: {m}\n'.format(s=SCRIPT_NAME, m=msg))

def outcome(target, return_code, dropped_inc, dropped_exc):
    return {'target': target, 
            'reason_code': return_code,
            'dropped_inc': dropped_inc,
            'dropped_exc': dropped_exc}
def taxa_in_tree(tree, taxa_list):
    t = []
    dropped = []
    for i in taxa_list:
        ind = tree.label2index.get(i)
        if ind is not None:
            t.append(tree.taxon_namespace[ind])
        else:
            dropped.append(i)
    return t, dropped
def bits_in_tree(tree, taxa_list):
    t, dropped = taxa_in_tree(tree, taxa_list)
    return [tree.label2bit[i.label] for i in t], dropped
def add_node_based_phyloreferenced_annotation(tree, annotation):
    in_tree, dropped_inc = taxa_in_tree(tree, annotation.des)
    if len(in_tree) == 0:
        return outcome(False, Reason.NO_INC_DESIGNATORS_IN_TREE, dropped_inc, None)
    if len(in_tree) == 1:
        mrca = tree.find_node_with_taxon_label(in_tree[0].label)
    else:
        mrca = tree.mrca(taxa=in_tree)
    assert mrca is not None
    return outcome(mrca, Reason.SUCCESS, None, None)

def add_stem_based_phylorefenced_annotation(tree, annotation):
    in_tree, dropped_inc = taxa_in_tree(tree, annotation.des)
    exclude, dropped_ex = bits_in_tree(tree, annotation.exclude_ancs_of)
    if not exclude:
        return outcome(tree.seed_node, Reason.SUCCESS, dropped_inc, dropped_ex)
    exc_bit_set = 0
    for e in exclude:
        exc_bit_set |= e
    if len(in_tree) == 0:
        return outcome(False, Reason.NO_INC_DESIGNATORS_IN_TREE, dropped_inc, dropped_ex)
    if len(in_tree) == 1:
        mrca = tree.find_node_with_taxon_label(in_tree[0].label)
    else:
        mrca = tree.mrca(taxa=in_tree)
    assert mrca is not None
    if mrca.edge.split_bitmask & exc_bit_set:
        return outcome(False, Reason.MRCA_HAS_EXCLUDED, dropped_inc, dropped_ex)
    deepest_valid = mrca
    curr = mrca.parent_node
    while curr is not None and not (curr.edge.split_bitmask & exc_bit_set):
        deepest_valid = curr
        curr = curr.parent_node
    return outcome(curr.edge, Reason.SUCCESS, dropped_inc, dropped_ex)

def add_phyloreferenced_annotation(tree, annotation):
    if annotation.exclude_ancs_of:
        r = add_stem_based_phylorefenced_annotation(tree, annotation)
    else:
        r = add_node_based_phyloreferenced_annotation(tree, annotation)
    if r['target']:
        r['target'].phylo_ref.append(annotation)
        annotation.applied_to.append((tree, r['target']))
    return r


class PhyloReferencedAnnotation(object):
    def __init__(self, serialized):
        self.des = [str(i) for i in serialized['descendants']]
        self.exclude_ancs_of = [str(i) for i in serialized.get('excludes_ancestors_of', [])]
        self.applied_to = []
if __name__ == '__main__':
    try:
        tree_file, annotations_file = sys.argv[1:]
    except:
        sys.exit('expecting 2 arguments: a tree file and a JSON file of annotations')
    tree_list = dendropy.TreeList.get_from_path(tree_file,
                                                'newick',
                                                suppress_internal_node_taxa=False)
    for tree in tree_list:
        tree.print_plot()
        mod_encode_splits(tree, delete_outdegree_one=False, internal_node_taxa=True)
        tree.label2index = {}
        tree.label2bit = {}
        curr_bit = 1
        for n, taxon in enumerate(tree.taxon_namespace):
            assert taxon.label not in tree.label2index
            tree.label2index[taxon.label] = n
            tree.label2bit[taxon.label] = curr_bit
            curr_bit <<= 1
        print tree.label2index
        print tree.label2bit
        for node in tree.preorder_node_iter():
            node.phylo_ref = []
            if node.edge:
                node.edge.phylo_ref = []
            print node.edge.split_bitmask
            if node.taxon:
                print node.taxon.label

    a_f = codecs.open(annotations_file, 'rU', encoding='utf-8')
    annot_list = json.load(a_f)
    if not isinstance(annot_list, list):
        annot_list = [annot_list]

    for tree_index, tree in enumerate(tree_list):
        num_tried = 0
        num_added = 0
        for annot_index, annotation in enumerate(annot_list):
            a = PhyloReferencedAnnotation(annotation)
            response = add_phyloreferenced_annotation(tree, a)
            x = response['target']
            if x:
                num_added += 1
            else:
                debug('Annotation {a} could not be added to tree {t}'.format(a=annot_index, t=tree_index))
            num_tried += 1
        debug('{a}/{t} annotations added to tree {i}'.format(a=num_added, t=num_tried, i=tree_index))
