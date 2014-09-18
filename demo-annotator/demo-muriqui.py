#!/usr/bin/env python
import dendropy
from dendropy.utility import container
from peyotl.api import APIWrapper
import codecs
import json
import sys
import os
TAXOMACHINE = APIWrapper().taxomachine
TREEMACHINE = APIWrapper().tree_of_life
SCRIPT_NAME = os.path.split(sys.argv[0])[1]
class Reason(object):
    NO_INC_DESIGNATORS_IN_TREE = 0
    SUCCESS = 1
    MRCA_HAS_EXCLUDED = 2
    ERROR_CHECK_FAILED = 3
    def to_str(c):
        if c == Reason.ERROR_CHECK_FAILED:
            return 'An error check failed'
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

class MappingOutcome(object):
    def __init__(self, attached_to=None, reason_code=Reason.SUCCESS, missing_inc=tuple(), missing_exc=tuple()):
        self.attached_to = attached_to
        self.reason_code = reason_code
        self.missing_inc = missing_inc
        self.missing_exc = missing_exc
        self.failed_error_checks = []
        self.failed_warning_checks = []
    def add_failed_error_check(self, check):
        self.failed_error_checks.append(check)
        self.reason_code = Reason.ERROR_CHECK_FAILED
    def add_failed_warning_check(self, check):
        self.failed_warning_checks.append(check)
    def explain(self):
        if self.reason_code == Reason.SUCCESS:
            return 'succcessfully mapped to {}'.format(self.attached_to)
        if self.reason_code == Reason.ERROR_CHECK_FAILED:
            return 'Error check ({}) failed.'.format(self.failed_error_checks[0].explain())
        return 'Attaching the annotation to the tree failed ({})'.format(Reason.to_str(self.reason_code))



def taxa_in_tree(tree, taxa_list, bits=False):
    t = []
    dropped = []
    for i in taxa_list:
        ind = tree.label2index.get(i)
        if ind is not None:
            if bits:
                t.append(tree.label2bit[i])
            else:
                t.append(tree.taxon_namespace[ind])
        else:
            dropped.append(i)
    return t, dropped
def bits_in_tree(tree, taxa_list):
    t, dropped = taxa_in_tree(tree, taxa_list)
    return [tree.label2bit[i.label] for i in t], dropped
def find_node_based_target(tree, annotation):
    resp = _find_mrca_and_verify_monophyly(tree, annotation)
    if isinstance(resp, MappingOutcome):
        return resp
    mrca, exc_bit_set, dropped_inc, dropped_exc = resp
    return MappingOutcome(mrca, Reason.SUCCESS, None, None)

def _find_mrca_and_verify_monophyly(tree, annotation):
    in_tree, dropped_inc = taxa_in_tree(tree, annotation.des)
    exclude, dropped_exc = bits_in_tree(tree, annotation.exclude_ancs_of)
    if not exclude:
        return MappingOutcome(tree.seed_node, Reason.SUCCESS, dropped_inc, dropped_exc)
    exc_bit_set = 0
    for e in exclude:
        exc_bit_set |= e
    if len(in_tree) == 0:
        return MappingOutcome(False, Reason.NO_INC_DESIGNATORS_IN_TREE, dropped_inc, dropped_exc)
    if len(in_tree) == 1:
        mrca = tree.find_node_with_taxon_label(in_tree[0].label)
    else:
        mrca = tree.mrca(taxa=in_tree)
    assert mrca is not None
    if mrca.edge.split_bitmask & exc_bit_set:
        return MappingOutcome(False, Reason.MRCA_HAS_EXCLUDED, dropped_inc, dropped_exc)
    return mrca, exc_bit_set, dropped_inc, dropped_exc

def find_stem_based_phylorefenced_annotation(tree, annotation):
    resp = _find_mrca_and_verify_monophyly(tree, annotation)
    if isinstance(resp, MappingOutcome):
        return resp
    mrca, exc_bit_set, dropped_inc, dropped_exc = resp
    deepest_valid = mrca
    curr = mrca.parent_node
    while (curr is not None) and ((curr.edge.split_bitmask & exc_bit_set) == 0):
        #print 'no intersection of', curr.edge.split_bitmask, exc_bit_set
        deepest_valid = curr
        curr = curr.parent_node
    #if curr:
    #    print 'intersection of', curr.edge.split_bitmask, exc_bit_set
    return MappingOutcome(deepest_valid.edge, Reason.SUCCESS, dropped_inc, dropped_exc)
class CheckOutcome(object):
    def __init__(self, passed, check):
        self.passed = passed
        self.check = check
def check_passes_result(check):
    return CheckOutcome(True, check)
def perform_check(tree, node_or_edge, check):
    if check.passes(tree, node_or_edge):
        return check_passes_result(check)
    return CheckOutcome(False, check)
def add_phyloreferenced_annotation(tree, annotation):
    #debug('Trying annotation {}'.format(annotation.annot_id))
    if annotation.rooted_by == GroupType.BRANCH:
        r = find_stem_based_phylorefenced_annotation(tree, annotation)
    else:
        assert annotation.rooted_by == GroupType.NODE
        r = find_node_based_target(tree, annotation)
    if r.reason_code != Reason.SUCCESS:
        return r
    for check in annotation.error_checks:
        check_result = perform_check(tree, r.attached_to, check)
        if not check_result.passed:
            r.add_failed_error_check(check)
            return r
    for check in annotation.warning_checks:
        check_result = perform_check(tree, r.attached_to, check)
        if not check_result.passed:
            r.add_failed_warning_check(check)
    #debug('Adding annotation {a} to {e}'.format(a=annotation.annot_id, e=r.attached_to))
    r.attached_to.phylo_ref.append(annotation)
    annotation.applied_to.append((tree, r.attached_to))
    assert r.reason_code == Reason.SUCCESS
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
class _CheckBase(object):
    pass
def expand_clade_using_ott(ott_id):
    # TAXOMACHINE.subtree(ott_id)
    return [ott_id]
class MonophylyCheck(object):
    def __init__(self, *valist):
        self.clade_list = [str(i) for i in valist]
    def explain(self):
        return 'REQUIRE_MONOPHYLETIC({})'.format(', '.join(self.clade_list))
    def passes(self, tree, node_or_edge):
        bitmask = 0
        for c in self.clade_list:
            expanded = expand_clade_using_ott(c)
            in_tree, m = taxa_in_tree(tree, expanded, bits=True)
            for x in in_tree:
                bitmask |= x
        if (bitmask == 0) or (bitmask not in tree.split_edges):
            return False
        return True
class CladeExcludesCheck(object):
    def __init__(self, *valist):
        self.clade_list = [str(i) for i in valist]
        self.failed = None
    def explain(self):
        return 'TARGET_EXCLUDES({})'.format(', '.join(self.clade_list))
    def passes(self, tree, node_or_edge):
        self.failed = None
        if not node_or_edge:
            return None
        try:
            edge = node_or_edge.edge
        except:
            edge = node_or_edge
        for c in self.clade_list:
            expanded = expand_clade_using_ott(c)
            in_tree, m = taxa_in_tree(tree, expanded, bits=True)
            exc_code = 0
            for t in in_tree:
                exc_code |= t
            if exc_code * edge.split_bitmask:
                self.failed = c
                return False
        return True

_CHECK_CODE_TO_TYPE = {'REQUIRE_MONOPHYLETIC': MonophylyCheck, 
                       'TARGET_EXCLUDES': CladeExcludesCheck, }
def deserialize_check(from_json):
    type_code = from_json[0]
    t = _CHECK_CODE_TO_TYPE[type_code]
    return t(*from_json[1:])

class PhyloReferencedAnnotation(object):
    def __init__(self, serialized):
        self.annot_id = serialized['_id'] # unique ID
        self.target = serialized['oa:hasTarget']
        self.annotated_at = serialized['oa:annotatedAt']
        self.annotated_by = serialized['oa:annotatedBy']
        self.body = serialized['oa:hasBody']
        self.des = [str(i) for i in self.target['included_ids']]
        self.exclude_ancs_of = [str(i) for i in self.target.get('excluded_ids', [])]
        self.rooted_by = GroupType.to_code(self.target['type'])
        self.error_checks = [deserialize_check(i) for i in self.target.get('error_checks', [])]
        self.warning_checks = [deserialize_check(i) for i in self.target.get('warning_checks', [])]

        self.applied_to = []
    def get_summary(self):
        return json.dumps(self.serialize())
    summary = property(get_summary)
    def serialize(self):
        return {'oa:hasTarget': self.target,
                'oa:annotatedBy': self.annotated_by,
                'oa:annotatedAt': self.annotated_at,
                'oa:hasBody': self.body,}
def all_numeric_taxa(tree_list):
    for tree in tree_list:
        for taxon in tree.taxon_namespace:
            try:
                int(taxon.label)
            except:
                return False
    return True
def convert_taxon_labels_to_ott_id(tree_list):
    tree  = tree_list[0]
    for taxon in tree.taxon_namespace:
        try:
            int(taxon.label)
        except:
            s = taxon.label.split('_')
            try:
                if len(s) < 2:
                    s = taxon.label.split(' ')
                if len(s) < 2:
                    sys.stderr.write('name without underscore "{}" found.\n'.format(taxon.label))
                    assert False

                o = s[-1]
                if not o.startswith('ott'):
                    o = taxon.label.split(' ')[-1]
                    if not o.startswith('ott'):
                        sys.stderr.write('name without trailing _ott# "{}" found.\n'.format(taxon.label))
                        assert False
                ott_id = o[3:]
                taxon.label = ott_id
            except:
                sys.exit('Currently the tree must be either labelled with only ott IDs or the using the name_ott<OTTID> convention.')
UNNAMED_NODE_COUNT = 0
def get_node_out_id(node):
    global UNNAMED_NODE_COUNT
    try:
        if node.label:
            return node.label
    except:
        pass
    if node.taxon and node.taxon.label:
        return node.taxon.label
    l = 'AUTOGENID' + str(UNNAMED_NODE_COUNT)
    UNNAMED_NODE_COUNT += 1
    node.label = l
    return l
def main(tree_filename, annotations_filename, out_tree_file_obj, out_table_file_obj):
    if not os.path.exists(tree_filename):
        raise ValueError('tree file "{}" does not exist'.format(tree_filename))
    tree_list = dendropy.TreeList.get_from_path(tree_filename,
                                                'newick',
                                                suppress_internal_node_taxa=False)
    if len(tree_list) != 1:
        sys.exit('sorry, the tree input can only be one tree, now')
    if not all_numeric_taxa(tree_list):
        convert_taxon_labels_to_ott_id(tree_list)
    for tree in tree_list:
        #tree.print_plot(show_internal_node_ids=True)
        mod_encode_splits(tree, delete_outdegree_one=False, internal_node_taxa=True)
        tree.label2index = {}
        tree.label2bit = {}
        curr_bit = 1
        for n, taxon in enumerate(tree.taxon_namespace):
            assert taxon.label not in tree.label2index
            tree.label2index[taxon.label] = n
            tree.label2bit[taxon.label] = curr_bit
            curr_bit <<= 1
        #print tree.label2index
        #print tree.label2bit
        for node in tree.preorder_node_iter():
            node.phylo_ref = []
            if node.edge:
                node.edge.phylo_ref = []
            #print node.edge.split_bitmask
            #if node.taxon:
            #   print node.taxon.label

    a_f = codecs.open(annotations_filename, 'rU', encoding='utf-8')
    annot_list = json.load(a_f)
    if not isinstance(annot_list, list):
        annot_list = [annot_list]
    for tree_index, tree in enumerate(tree_list):
        num_tried = 0
        num_added = 0
        unadded = []
        for annotation in annot_list:
            a = PhyloReferencedAnnotation(annotation)
            response = add_phyloreferenced_annotation(tree, a)
            if response.reason_code == Reason.SUCCESS:
                num_added += 1
            else:
                unadded.append((a, response))
                #debug('Annotation {a} could not be added to tree {t}'.format(a=a.annot_id, t=tree_index))
            num_tried += 1
        debug('{a}/{t} annotations added to tree {i}'.format(a=num_added, t=num_tried, i=tree_index))
        # Report tree and annotations
        tree.write(out_tree_file_obj, 'newick', node_label_compose_func=get_node_out_id)
        out_table_file_obj.write('type\taddress\tannot-id\n')
        for node in tree.preorder_node_iter():
            if node.taxon:
                l = node.taxon.label
            else:
                l = '@' + str(id(node))
            if node.phylo_ref:
                for a in node.phylo_ref:
                    out_table_file_obj.write('node\t{n}\t{a}\n'.format(n=get_node_out_id(node), a=a.annot_id))
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
                if e.phylo_ref:
                    for a in e.phylo_ref:
                        out_table_file_obj.write('edge\t{n}\t{a}\n'.format(n=get_node_out_id(node), a=a.annot_id))
        # Report unadded annotations
        for annotation, add_record in unadded:
            out_table_file_obj.write('NA\t\t{a}\n'.format(a=a.annot_id))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('demo of muriqui mapping annotations to a tree')
    parser.add_argument('--taxon-tree', type=int, help='ott ID that will be used to fetch taxonomy/subtree to use as the tree to be annotated')
    parser.add_argument('--tree-node', type=int, help='treemachine node ID that will be used to fetch tree_of_life/subtree to use as the tree to be annotated')
    parser.add_argument('--tree-ott', type=int, help='ott ID that will be used to fetch tree_of_life/subtree to use as the tree to be annotated')
    parser.add_argument('--tree-file', help='filepath to newick file with labels as ott IDs or using the name_ott#### convention')
    parser.add_argument('--out-table', required=True, help='file to output with the annotation placements')
    parser.add_argument('--out-tree', required=True, help='file to output with a tree with IDs to be used with the out-table')
    parser.add_argument('json', help='filepath to JSON file with annotations')
    args = parser.parse_args()
    annotations_file = args.json
    o_tree = open(args.out_tree, 'w')
    o_table = open(args.out_table, 'w')

    if args.tree_file is not None:
        tree_file = args.tree_file
    else:
        if args.taxon_tree is not None:
            resp = TAXOMACHINE.subtree(int(args.taxon_tree))['subtree']
        elif args.tree_node is not None:
            resp = TREEMACHINE.subtree(node_id=int(args.tree_node))['newick']
        elif args.tree_ott is not None:
            resp = TREEMACHINE.subtree(ott_id=int(args.tree_ott))['newick']
        else:
            sys.exit('must specify a tree\n')
        import tempfile
        handle, tmpf = tempfile.mkstemp(suffix='.tre')
        sys.stderr.write('Writing fetched tree to "{}"'.format(tmpf))
        os.close(handle)
        o = codecs.open(tmpf, 'w', encoding='utf-8')
        o.write(resp)
        o.write(';\n')
        o.close()
        tree_file = tmpf
    main(tree_file, annotations_file, o_tree, o_table)
