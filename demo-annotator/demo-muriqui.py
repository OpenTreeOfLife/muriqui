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
    def to_str(c):
        if c == Reason.NO_INC_DESIGNATORS_IN_TREE:
            return 'no specifiers to be included were in the tree'
        if c == Reason.MRCA_HAS_EXCLUDED:
            return 'the include group is paraphyletic with respect to member/members of the exclude group'
        return ''
    to_str = staticmethod(to_str)

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
                assert p
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
    resp = _find_mrca_and_verify_monophyly(tree, annotation)
    if isinstance(resp, dict):
        return resp
    mrca, exc_bit_set, dropped_inc, dropped_exc = resp
    return outcome(mrca, Reason.SUCCESS, None, None)

def _find_mrca_and_verify_monophyly(tree, annotation):
    in_tree, dropped_inc = taxa_in_tree(tree, annotation.des)
    exclude, dropped_exc = bits_in_tree(tree, annotation.exclude_ancs_of)
    if not exclude:
        return outcome(tree.seed_node, Reason.SUCCESS, dropped_inc, dropped_exc)
    exc_bit_set = 0
    for e in exclude:
        exc_bit_set |= e
    if len(in_tree) == 0:
        return outcome(False, Reason.NO_INC_DESIGNATORS_IN_TREE, dropped_inc, dropped_exc)
    if len(in_tree) == 1:
        mrca = tree.find_node_with_taxon_label(in_tree[0].label)
    else:
        mrca = tree.mrca(taxa=in_tree)
    assert mrca is not None
    if mrca.edge.split_bitmask & exc_bit_set:
        return outcome(False, Reason.MRCA_HAS_EXCLUDED, dropped_inc, dropped_exc)
    return mrca, exc_bit_set, dropped_inc, dropped_exc
def add_stem_based_phylorefenced_annotation(tree, annotation):
    resp = _find_mrca_and_verify_monophyly(tree, annotation)
    if isinstance(resp, dict):
        return resp
    mrca, exc_bit_set, dropped_inc, dropped_exc = resp
    deepest_valid = mrca
    curr = mrca.parent_node
    while (curr is not None) and ((curr.edge.split_bitmask & exc_bit_set) == 0):
        print 'no intersection of', curr.edge.split_bitmask, exc_bit_set
        deepest_valid = curr
        curr = curr.parent_node
    if curr:
        print 'intersection of', curr.edge.split_bitmask, exc_bit_set
    return outcome(deepest_valid.edge, Reason.SUCCESS, dropped_inc, dropped_exc)

def add_phyloreferenced_annotation(tree, annotation):
    if annotation.rooted_by == GroupType.BRANCH:
        r = add_stem_based_phylorefenced_annotation(tree, annotation)
    else:
        assert annotation.rooted_by == GroupType.NODE
        r = add_node_based_phyloreferenced_annotation(tree, annotation)
    if r['target']:
        r['target'].phylo_ref.append(annotation)
        annotation.applied_to.append((tree, r['target']))
    return r

class GroupType:
    BRANCH, NODE = range(2)
    def to_str(c):
        if c == GroupType.BRANCH:
            return 'branch'
        assert c == GroupType.NODE
        return 'node'
    to_str = staticmethod(to_str)
    def to_code(c):
        if c.lower() == 'branch':
            return GroupType.BRANCH
        assert c.lower() == 'node'
        return GroupType.NODE
    to_code = staticmethod(to_code)

class PhyloReferencedAnnotation(object):
    def __init__(self, serialized):
        self.target = serialized['oa:hasTarget']
        self.annotated_at = serialized['oa:annotatedAt']
        self.annotated_by = serialized['oa:annotatedBy']
        self.body = serialized['oa:hasBody']
        self.des = [str(i) for i in self.target['included_ids']]
        self.exclude_ancs_of = [str(i) for i in self.target.get('excluded_ids', [])]
        self.rooted_by = GroupType.to_code(self.target['type'])
        self.applied_to = []
    def get_summary(self):
        return json.dumps(self.serialize())
    summary = property(get_summary)
    def serialize(self):
        return {'oa:hasTarget': self.target,
                'oa:annotatedBy': self.annotated_by,
                'oa:annotatedAt': self.annotated_at,
                'oa:hasBody': self.body,}
def main(tree_filename, annotations_filename):
    tree_list = dendropy.TreeList.get_from_path(tree_filename,
                                                'newick',
                                                suppress_internal_node_taxa=False)
    for tree in tree_list:
        tree.print_plot(show_internal_node_ids=True)
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

    a_f = codecs.open(annotations_filename, 'rU', encoding='utf-8')
    annot_list = json.load(a_f)
    if not isinstance(annot_list, list):
        annot_list = [annot_list]

    for tree_index, tree in enumerate(tree_list):
        num_tried = 0
        num_added = 0
        unadded = []
        for annot_index, annotation in enumerate(annot_list):
            a = PhyloReferencedAnnotation(annotation)
            response = add_phyloreferenced_annotation(tree, a)
            x = response['target']
            if x:
                num_added += 1
            else:
                unadded.append((a, response))
                debug('Annotation {a} could not be added to tree {t}'.format(a=annot_index, t=tree_index))
            num_tried += 1
        debug('{a}/{t} annotations added to tree {i}'.format(a=num_added, t=num_tried, i=tree_index))
        # Report tree and annotations
        tree.print_plot(show_internal_node_ids=True)
        for node in tree.preorder_node_iter():
            if node.taxon:
                l = node.taxon.label
            else:
                l = '@' + str(id(node))
            print 'Node "{l}":'.format(l=l)
            if node.phylo_ref:
                for a in node.phylo_ref:
                    print '   ', a.summary
            else:
                print '    <NO ANNOTATIONS>'
            e = node.edge
            if e:
                p = e.tail_node
                if p:
                    if p.taxon:
                        pl = p.taxon.label
                    else:
                        pl = '@' + str(id(p))
                else:
                    pl = 'None'
                print 'Edge "{p}" -> "{l}"'.format(p=pl, l=l)
                if e.phylo_ref:
                    for a in e.phylo_ref:
                        print '   ', a.summary
                else:
                    print '    <NO ANNOTATIONS>'
        # Report unadded annotations
        if len(unadded) > 0:
            print 'Unattached annotations'
        for annotation, add_record in unadded:
            print 'reason={r}. annotation={a}'.format(a=annotation.summary,
                                                      r=Reason.to_str(add_record['reason_code']))


if __name__ == '__main__':
    try:
        tree_file, annotations_file = sys.argv[1:]
    except:
        sys.exit('expecting 2 arguments: a tree file and a JSON file of annotations')
    main(tree_file, annotations_file)
