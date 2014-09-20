#!/usr/bin/env python
import dendropy
from dendropy.utility import container
from peyotl.api import APIWrapper
from cStringIO import StringIO
import codecs
import datetime
import json
import os
import random
import shutil
import string
import sys
import time
import unittest
TAXOMACHINE = APIWrapper().taxomachine
TREEMACHINE = APIWrapper().tree_of_life
SCRIPT_NAME = os.path.split(sys.argv[0])[1]
class Reason(object):
    NO_INC_DESIGNATORS_IN_TREE = 0
    SUCCESS = 1
    MRCA_HAS_EXCLUDED = 2
    ERROR_CHECK_FAILED = 3
    def to_str(c):
        if c == Reason.SUCCESS:
            return 'success'
        if c == Reason.ERROR_CHECK_FAILED:
            return 'an error check failed'
        if c == Reason.NO_INC_DESIGNATORS_IN_TREE:
            return 'no specifiers to be included were in the tree'
        if c == Reason.MRCA_HAS_EXCLUDED:
            return 'the include group is paraphyletic with respect to member/members of the exclude group'
        return ''
    to_str = staticmethod(to_str)

def debug(msg):
    sys.stderr.write('{s}: {m}\n'.format(s=SCRIPT_NAME, m=msg))

class MappingOutcome(object):
    def __init__(self, attached_to, reason_code, missing_inc, missing_exc):
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

class TargetTree(object):

    tree = None
    _unnamed_node_count = 0

    unadded_annotations = []
    _num_tried = 0
    _num_added = 0

    _name_converter = None
    
    @property
    def number_annotations_tried(self):
        return self._num_tried
    @property
    def number_annotations_added(self):
        return self._num_added
    
    def __init__(self, tree, use_taxonomy=True):
        self.tree = tree
        self.preorder_node_iter = self.tree.preorder_node_iter
        self.print_plot = self.tree.print_plot
        self.write = self.tree.write
        self._use_taxonomy = use_taxonomy

        if use_taxonomy:
            self._name_converter = OTTNameConverter()
            # check if all taxa have numeric names (assumed to be ott ids)
            all_numeric_taxa = True
            for taxon in self.tree.taxon_namespace:
                try:
                    int(taxon.label)
                except:
                    all_numeric_taxa = False
                    break

            # assuming we're dealing with ott ids here
            if not all_numeric_taxa:
                for taxon in self.tree.taxon_namespace:
                    try:
                        int(taxon.label)
                    except:
                        ott_id = self._name_converter.concat_taxon_label_to_ott_id(taxon.label)
                        taxon.label = ott_id

        #tree.print_plot(show_internal_node_ids=True)
        self.mod_encode_splits(tree, delete_outdegree_one=False, internal_node_taxa=True)
        self.tree.label2index = {}
        self.tree.label2bit = {}
        curr_bit = 1
        for n, taxon in enumerate(tree.taxon_namespace):
            assert taxon.label not in tree.label2index
            self.tree.label2index[taxon.label] = n
            self.tree.label2bit[taxon.label] = curr_bit
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
            
    def get_taxa_in_tree(self, ids, bits=False):
        if isinstance(ids,str):
            ids = [ids,]
        if self._use_taxonomy:
            ids = self._expand_ids(ids)
        found = []
        not_found = []
        for i in ids:
#            debug("searching for {} in tree".format(i))
            ind = self.tree.label2index.get(i)
            if ind is not None:
                if bits:
                    found.append(self.tree.label2bit[i])
                else:
                    found.append(self.tree.taxon_namespace[ind])
            else:
                not_found.append(i)
        return found, not_found

    def get_bits_in_tree(self, ids):
        if isinstance(ids,str):
            ids = [ids,]
        if self._use_taxonomy:
            ids = self._expand_ids(ids)
        found, not_found = self.get_taxa_in_tree(ids)
        return [self.tree.label2bit[i.label] for i in found], not_found

    def get_mrca(self, taxa):
        mrca = None
        if len(taxa) == 1:
            mrca = self.tree.find_node_with_taxon_label(taxa[0].label)
        else:
            mrca = self.tree.mrca(taxa=taxa)
        return mrca

    def _expand_ids(self, ids):
        e = []
        for i in ids:
            e.extend(self._name_converter.expand_clade_using_ott(i))
        return e
            
    def find_node_based_target(self, annotation):

        # get included nodes. if none in tree, fail
        included, dropped_inc = self.get_taxa_in_tree(annotation.target.ids_to_include)
        if len(included) < 1:
            return MappingOutcome(None, Reason.NO_INC_DESIGNATORS_IN_TREE, dropped_inc, None)

        mrca = self.get_mrca(included)
        assert mrca != None

        return MappingOutcome(mrca, Reason.SUCCESS, dropped_inc, None)

    def find_stem_based_target(self, annotation):

        # get included nodes. if none in tree, fail
        included, dropped_inc = self.get_taxa_in_tree(annotation.target.ids_to_include)
        if len(included) < 1:
            return MappingOutcome(None, Reason.NO_INC_DESIGNATORS_IN_TREE, dropped_inc, None)

        # get excluded bits. if none in tree, return root
        excluded, dropped_exc = self.get_bits_in_tree(annotation.target.ids_to_exclude)
        if len(excluded) < 1:
            return MappingOutcome(self.tree.seed_node, Reason.SUCCESS, dropped_inc, dropped_exc)

        # tally the excluded tips in bitset
        exc_bit_set = 0
        for e in excluded:
            exc_bit_set |= e

        # fail if the intersection of excluded bitset and mrca is not empty
        if (mrca.edge.split_bitmask & exc_bit_set) != 0:
            return MappingOutcome(None, Reason.MRCA_HAS_EXCLUDED, dropped_inc, dropped_exc)

        # get the mrca of the included nodes
        mrca = self.get_mrca(included)
        assert mrca is not None

        # find the deepest valid mrca that doesn't include any excluded nodes
        deepest_valid = mrca
        curr = mrca.parent_node
        while (curr is not None) and ((curr.edge.split_bitmask & exc_bit_set) == 0):
            #print 'no intersection of', curr.edge.split_bitmask, exc_bit_set
            deepest_valid = curr
            curr = curr.parent_node
        #if curr:
        #    print 'intersection of', curr.edge.split_bitmask, exc_bit_set
        return MappingOutcome(deepest_valid.edge, Reason.SUCCESS, dropped_inc, dropped_exc)

    def mod_encode_splits(self, create_dict=True, delete_outdegree_one=True, internal_node_taxa=False):
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
        taxon_namespace = self.tree.taxon_namespace
        if taxon_namespace is None:
            taxon_namespace = self.tree.infer_taxa()
        if create_dict:
            self.split_edges = {}
            split_map = self.split_edges
            # if tree.is_rooted:
            #     tree.split_edges = {}
            # else:
            #     atb = taxon_namespace.all_taxa_bitmask()
            #     d = container.NormalizedBitmaskDict(mask=atb)
            #     tree.split_edges = d
            # split_map = tree.split_edges
        if not self.tree.seed_node:
            return

        if delete_outdegree_one:
            sn = self.tree.seed_node
            if not self.tree.is_rooted:
                if len(sn._child_nodes) == 2:
                    self.tree.deroot()
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

        for edge in self.tree.postorder_edge_iter():
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
        if not self.tree.is_rooted:
            mask = self.tree.seed_node.edge.split_bitmask
            d = container.NormalizedBitmaskDict(mask=mask)
            for k, v in self.split_edges.items():
                d[k] = v
            self.split_edges = d

    def perform_check(self, node_or_edge, check):
        if check.passes(self, node_or_edge):
            return CheckOutcome(True, check)
        return CheckOutcome(False, check)

    def add_phyloreferenced_annotation(self, annotation):
        self._num_tried += 1
        
        # find the target
        if annotation.target.type == TargetType.BRANCH:
            r = self.find_stem_based_target(annotation)
        else:
            assert annotation.target.type == TargetType.NODE
            r = self.find_node_based_target(annotation)

        # if no target was found, fail
        if r.reason_code != Reason.SUCCESS:
            self.unadded_annotations.append((annotation, r))
            debug(self.unadded_annotations)
            return r

        # check error conditions
        for check in annotation.target.error_checks:
            check_result = self.perform_check(r.attached_to, check)
            if not check_result.passed:
                r.add_failed_error_check(check)
                
                # stop if we fail an error check
                self.unadded_annotations.append((annotation, r))
                return r

        # check warning conditions
        for check in annotation.target.warning_checks:
            check_result = self.perform_check(r.attached_to, check)
            if not check_result.passed:
                r.add_failed_warning_check(check)

        # add the annotation
        r.attached_to.phylo_ref.append(annotation)
        annotation.applied_to.append((self.tree, r.attached_to))
        self._num_added += 1
        return r

    def get_node_out_id(self,node):
        try:
            if node.label:
                return node.label
        except:
            pass
        if node.taxon and node.taxon.label:
            return node.taxon.label
        l = 'AUTOGENID' + str(self._unnamed_node_count)
        self._unnamed_node_count += 1
        node.label = l
        return l

    def write_labeled_tree(self, tree_file):
        self.tree.write(tree_file, 'newick', node_label_compose_func=self.get_node_out_id)

    def write_table(self, table_file):
        table_file.write('type\ttarget_id\tannotation_id\treason\n')
        for node in self.tree.preorder_node_iter():
            if node.phylo_ref:
                for a in node.phylo_ref:
                    table_file.write('node\t{n}\t{a}\t{o}\n'.format(n=self.get_node_out_id(node), \
                            a=a.id, o=Reason.to_str(Reason.SUCCESS)))
            e = node.edge
            if e:
                if e.phylo_ref:
                    for a in e.phylo_ref:
                        table_file.write('edge\t{n}\t{a}\t{o}\n'.format(n=self.get_node_out_id(node), \
                                a=a.id, o=Reason.to_str(Reason.SUCCESS)))

        # report unadded annotations
        for annotation, result in self.unadded_annotations:
            table_file.write('NA\tNA\t{a}\t{o}\n'.format(a=annotation.id,\
                    o=Reason.to_str(result.reason_code)))

class CheckOutcome(object):
    def __init__(self, passed, check):
        self.passed = passed
        self.check = check

class TargetType(object):
    BRANCH, NODE = range(2)
    def to_str(c):
        if c == TargetType.BRANCH:
            return 'branch'
        assert c == TargetType.NODE
        return 'node'
    to_str = staticmethod(to_str)
    def to_code(c):
        if c.lower() == 'branch':
            return TargetType.BRANCH
        assert c.lower() == 'node'
        return TargetType.NODE
    to_code = staticmethod(to_code)

class OTTNameConverter(object):

    def __init__(self):
        self._EXP_CACHE = {}

    def get_ott_ids_from_taxon_namespace(self, ns):
        r = []
        for taxon in ns:
            x = self.concat_taxon_label_to_ott_id(taxon.label)
            r.append(x)
        return r

    def expand_clade_using_ott(self, ott_id):
        if ott_id in self._EXP_CACHE:
            return self._EXP_CACHE[ott_id]
        n = TAXOMACHINE.subtree(ott_id)['subtree']
        if n.startswith('('):
    #        n += ';'
            inp = StringIO(n)
            t = dendropy.Tree.get_from_stream(inp, 'newick')
            id_list = self.get_ott_ids_from_taxon_namespace(t.taxon_namespace)        
        else:
            id_list = [self.concat_taxon_label_to_ott_id(n)]
    
        # kludge to deal with dendropy bug where trailing semicolon is added to name
            id_list = [i.strip(";") for i in id_list]
    
        self._EXP_CACHE[ott_id] = id_list
        return id_list
    
    def concat_taxon_label_to_ott_id(self, label, from_taxom=False):
        s = label.split('_')
        try:
            if len(s) < 2:
                s = label.split(' ')
            if len(s) < 2:
                sys.stderr.write('name without underscore "{}" found.\n'.format(label))
                assert False
            o = s[-1]
            if not o.startswith('ott'):
                o = label.split(' ')[-1]
                if not o.startswith('ott'):
                    sys.stderr.write('name without trailing _ott# "{}" found.\n'.format(label))
                    assert False
            ott_id = o[3:]
        except:
            msg = 'Currently the tree must be either labelled with only ott IDs or the using the name_ott<OTTID> convention.'
            if from_taxom:
                msg += ' The tree was fetched from taxomachine internally to expand an ott ID.'
            raise ValueError(msg)
        return ott_id

    
class MonophylyCondition(object):
    def __init__(self, *valist):
        self.clade_list = [str(i) for i in valist]
    def explain(self):
        return 'REQUIRE_MONOPHYLETIC({})'.format(', '.join(self.clade_list))
    def passes(self, tree, node_or_edge):
        bitmask = 0
        for c in self.clade_list:
#            expanded = expand_clade_using_ott(c)
            in_tree = tree.get_taxa_in_tree(c, bits=True)[0]
#            in_tree = tree.taxa_in_tree(c, bits=True)[0]
            for x in in_tree:
                bitmask |= x
        if (bitmask == 0) or (bitmask not in tree.split_edges):
            return False
        return True
    def to_json(self):
        return ["REQUIRE_MONOPHYLETIC",] + self.clade_list

class TargetExcludesCondition(object):
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
#            expanded = expand_clade_using_ott(c)
            in_tree = tree.get_taxa_in_tree(c, bits=True)[0]
#            in_tree = tree.taxa_in_tree(c, bits=True)[0]
            exc_code = 0
            for t in in_tree:
                exc_code |= t
            if exc_code * edge.split_bitmask:
                self.failed = c
                return False
        return True
    def to_json(self):
        return ['TARGET_EXCLUDES',] + self.clade_list

class ReferenceCondition(object):

    _CODE_TO_TYPE = {
        'REQUIRE_MONOPHYLETIC': MonophylyCondition,
        'TARGET_EXCLUDES': TargetExcludesCondition,
    }

    @classmethod
    def from_data(cls,data):
        type_code = data[0]
        t = cls.get_type_from_code(type_code)
        return t(*data[1:])

    @staticmethod
    def get_type_from_code(code):
        return ReferenceCondition._CODE_TO_TYPE[code]
        
class ReferenceTarget(object):
    def __init__(self, group_type):
        self._type = TargetType.to_code(group_type)
        self._ids_to_include = []
        self._ids_to_exclude = []
        self._error_checks = []
        self._warning_checks = []
    @property
    def type(self):
        return self._type
    @type.setter
    def type(self, _type):
        self._type = _type
    @property
    def ids_to_include(self):
        return self._ids_to_include
    @property
    def ids_to_exclude(self):
        return self._ids_to_exclude
    @property
    def error_checks(self):
        return self._error_checks
    @property
    def warning_checks(self):
        return self._warning_checks

    def include_specifiers(self, specifiers):
        for s in specifiers:
            self._ids_to_include.append(s)
    def exclude_specifiers(self, specifiers):
        for s in specifiers:
            self._ids_to_exclude.append(s)

    def add_error_condition(self, condition):
        self._error_checks.append(condition)
    def add_warning_condition(self, condition):
        self._warning_checks.append(condition)
    @classmethod
    def from_data(cls, data):
        t = cls(data['type'])
        t.include_specifiers(data['included_ids'])
        t.exclude_specifiers(data.get('excluded_ids', []))
        for s in data.get('error_checks', []):
            t.add_error_condition(ReferenceCondition.from_data(s))
        for s in data.get('warning_checks', []):
            t.add_warning_condition(ReferenceCondition.from_data(s))
        return t
    def to_json(self):
        return {
            "type": TargetType.to_str(self._type),
            "included_ids": self._ids_to_include,
            "excluded_ids": self._ids_to_exclude,
            "error_checks": [x.to_json() for x in self._error_checks],
            "warning_checks": [x.to_json() for x in self._warning_checks],
        }
        
class Entity(object):
    def __init__(self):
        self._name = ""
        self._url = ""
        self._description = ""
        self._version = ""
        self._invocation = {}
        self._type = "prov:Entity"
    @property
    def url(self):
        return self._url
    @url.setter
    def url(self, url):
        self._url = url
    @property
    def description(self):
        return self._description
    @url.setter
    def description(self, description):
        self._description = description
    @property
    def version(self):
        return self._version
    @url.setter
    def version(self, version):
        self._version = version
    @property
    def invocation(self):
        return self._invocation
    @url.setter
    def invocation(self, invocation):
        self._invocation = invocation
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, name):
        self._name = name
    @property
    def type(self):
        return self._type
    @classmethod
    def from_data(cls, data):
        e = cls()
        e.name = data["name"]
        return e
    def to_json(self):
        return {
            "type": self._type,
            "name": self._name,
            "url": self._url,
            "description": self._description,
            "version": self._version,
            "invocation": self._invocation,
        }
        
class PhyloReferencedAnnotation(object):
    def __init__(self):
        self._id = None # unique ID
        self._target = None
        self._annotated_at = ""
        self._annotated_by = None
        self._body = None
        self._applied_to = []
    @property
    def id(self):
        return self._id
    @id.setter
    def id(self, _id):
        self._id = _id
    @property
    def target(self):
        return self._target
    @target.setter
    def target(self, target):
        self._target = target
    @property
    def annotated_at(self):
        return self._annotated_at
    @annotated_at.setter
    def annotated_at(self, datetime):
        self._annotated_at = datetime
    @property
    def annotated_by(self):
        return self._annotated_by
    @annotated_by.setter
    def annotated_by(self, entity):
        self._annotated_by = entity
    @property
    def body(self):
        return self._body
    @body.setter
    def body(self, body):
        self._body = body
    @property
    def applied_to(self):
        return self._applied_to
    @applied_to.setter
    def applied_to(self, applied_to):
        self._applied_to = applied_to
    @property
    def summary(self):
        return json.dumps(self.to_json())
    @classmethod
    def from_data(cls, data):
        a = cls()
        a.id = data['_id']
        a.target = ReferenceTarget.from_data(data['oa:hasTarget'])
        a.annotated_at = data['oa:annotatedAt']
        a.annotated_by = Entity.from_data(data['oa:annotatedBy'])
        a.body = data['oa:hasBody']
        return a
    def to_json(self):
        return {
            '_id': self._id,
            'oa:hasTarget': self.target.to_json(),
            'oa:annotatedBy': self.annotated_by.to_json(),
            'oa:annotatedAt': self.annotated_at,
            'oa:hasBody': self.body,
        }
        
class RandomAnnotation(PhyloReferencedAnnotation):
    
    MAX_FLOAT_VALUE = 100000.0
    MAX_INT_VALUE = 1000000000
    MAX_STRING_LENGTH = 40

    MAX_ITEMS = 5
    MAX_DEPTH = 5
    MAX_KEY_LENGTH = 40
    MAX_TARGET_LENGTH = 40
    MAX_ERROR_CONDITIONS = 5
    MAX_WARNING_CONDITIONS = 5

    simple_string_chars = string.letters + string.digits + string.punctuation

    def __init__(self, id, random_seed=None, use_utf8=False):
    
        if random_seed is not None:
            random.seed(random_seed)
        else:
            random.seed(time.clock()) 

        if use_utf8:
            self.get_random_string = self._get_random_string_utf8 
        else:
            self.get_random_string = self._get_random_string_ascii

        self.id = id

        e = Entity()
        e.name = self.get_random_string(random.randrange(30))
        e.url = "http://"+self.get_random_string(random.randrange(100))
        e.description = self.get_random_string(random.randrange(300))
        e.version = self.get_random_string(random.randrange(20))
        e.invocation = self._get_random_object(random.randrange(4))

        self.annotated_at = datetime.datetime.now().isoformat()
        self.annotated_by = e
        self.body = self._get_random_object(random.randrange(2,10))
            
        # create a random reference type
        self.target = ReferenceTarget(random.sample(["node","branch"],1)[0])

        # generate random included ids        
        included = set()
        n_to_include = random.randrange(1,self.MAX_TARGET_LENGTH)
        while len(included) < n_to_include:
            included.add(random.randrange(self.MAX_INT_VALUE))

        # generate random excluded ids (only for branch references)
        if self.target.type == TargetType.BRANCH:
            excluded = set()
            n_to_exclude = random.randrange(self.MAX_TARGET_LENGTH)
            while len(excluded) < n_to_exclude:
                t = random.randrange(self.MAX_INT_VALUE)
                if t not in included:
                    excluded.add(t)
            self.target.exclude_specifiers(excluded)

        self.target.include_specifiers(included)

        # add zero or more error checks
        for i in range(random.randrange(self.MAX_ERROR_CONDITIONS)):
            self.target.add_error_condition(self._get_random_condition(included))
    
        # add zero or more warning checks
        for i in range(random.randrange(self.MAX_WARNING_CONDITIONS)):
            self.target.add_warning_condition(self._get_random_condition(included))

    def _get_random_string_utf8(self,length):
        s = StringIO()
        i = 0
        while i < length:
            try:
                s.write(unicode(os.urandom(8), encoding='UTF-8'))
                i += 1
            except ValueError:
                continue
        return s.getvalue()

    def _get_random_string_ascii(self,length):
        s = StringIO()
        i = 0
        while i < length:
            s.write(random.sample(self.simple_string_chars,1)[0])
            i+=1
        return s.getvalue()

    def _get_random_float(self):
        return random.random() * random.randrange(self.MAX_FLOAT_VALUE)

    def _get_random_int(self):
        return random.randrange(self.MAX_INT_VALUE)

    def _get_random_value(self, depth=0):

        # halt deep recursion
        if depth > MAX_DEPTH:
            return self._get_random_primitive()

        depth += 1
        
        # generate containers infrequently
        r = random.randrange(8)
        if r == 0:
            value = []
            for i in range(random.randrange(self.MAX_ITEMS)):
                value.append(self._get_random_value(depth))
        elif r == 1:
            value = {}
            for i in range(random.randrange(self.MAX_ITEMS)):
                value[self.get_random_string(self.MAX_KEY_LENGTH)] = self._get_random_value(depth)
        else:
            # not a container, generate random primitive
            value = self._get_random_primitive()
    
        return value

    def _get_random_primitive(self):

        r = random.randrange(3)
        if r == 0:
            value = self.get_random_string(self.MAX_STRING_LENGTH)
        elif r == 1:
            value = self._get_random_float()
        else:
            value = self._get_random_int()

        return value

    def _get_random_object(self, number_of_elements):
        b = {}
        for i in range(number_of_elements):
            key = self.get_random_string(self.MAX_KEY_LENGTH)
            value = self._get_random_value()
            b[key] = value
        return b

    def _get_random_condition(self, included_specifiers):
        i = random.randrange(2)            
        # MonophylyCondition
        if i == 0: 
            specifiers = set()
            if len(included_specifiers) > 1: 
                n = random.randrange(1,len(included_specifiers))
            else:
                n = 1
            c = MonophylyCondition(*random.sample(included_specifiers,n))

        # TargetExcludesCondition
        elif i == 1:
            specifiers = set()
            n = random.randrange(1,self.MAX_TARGET_LENGTH)
            while len(specifiers) < n:
                r = random.randrange(self.MAX_INT_VALUE)
                if r not in included_specifiers:
                    specifiers.add(r)
            c = TargetExcludesCondition(*specifiers)
        return c

def main(tree_filename, annotations_filename, out_tree_file_path, out_table_file_path, use_taxonomy=True):
    
    # get the trees
    if not os.path.exists(tree_filename):
        raise ValueError('tree file "{}" does not exist'.format(tree_filename))
    tree_list = dendropy.TreeList.get_from_path(tree_filename,'newick',
            suppress_internal_node_taxa=False)
    if len(tree_list) < 1:
        sys.stderr.write('No trees in input list.')
        return False

    # get the annotations
    a_f = codecs.open(annotations_filename, 'rU', encoding='utf-8')
    annot_list = json.load(a_f)
    if not isinstance(annot_list, list):
        annot_list = [annot_list]
    annotations = []
    for a in annot_list:
        annotations.append(PhyloReferencedAnnotation.from_data(a))

    # annotate the trees
    for tree_index, t in enumerate(tree_list):
        tree = TargetTree(t, use_taxonomy=use_taxonomy)
        for a in annotations:
#            debug(a.summary)
            tree.add_phyloreferenced_annotation(a)

        # report tree and annotations
        with open(out_tree_file_path,"w") as out_tree_file: 
            tree.write_labeled_tree(out_tree_file)
        with open(out_table_file_path,"w") as out_table_file:
            tree.write_table(out_table_file)

class Tests(unittest.TestCase):

    def setUp(self):
        if not os.path.exists("tests"):
            os.mkdir("tests")
    
#    def tearDown(self):
#        shutil.rmtree("tests")
        
    def test_taxon_expansion(self):
        pass
        # test whether taxon id expansion is behaving itself
    
    def test_node_based_expansion(self):
        pass
    
    def test_branch_based_expansion(self):
        pass
    
    def test_canid_data(self):
        t="examples/canids.tre"
        a="examples/armadillo-annot.json" 
        ot="tests/canids-out-tree.tre","w"
        ob="tests/canids-out-table.tsv"
        
        main(t, a, ot, ob)
        
        out = []
        with open(ob) as outfile:
            for line in outfile:
                p = line.split()
                out.append((p[1],p[2]))
        out = tuple(out)
        self.failUnless(out == (("target_id","annotation_id"),("770319","3"),("770319","4"), \
                ("770319","5"),("770319","6"),("NA","1"),("NA","2")))
    
    def test_roundtrip_100_ascii_annotations_n_times(self):
        for i in range(100):
            self.roundtrip_random_annotation_n_times(random.randrange(1,10), False)

    def test_roundtrip_10_utf8_annotations_n_times(self):
        for i in range(10):
            self.roundtrip_random_annotation_n_times(random.randrange(1,10), True)
    
    def roundtrip_random_annotation_n_times(self, k, use_utf8=False):
        
            # make a random annotation, and store it as json
            r = RandomAnnotation(id=0, use_utf8=use_utf8)
            x = json.loads(r.summary)
            y = json.loads(r.summary)

            # roundtrip the json the specified number of times
            for i in range(k):
                try:
                    y = json.loads(PhyloReferencedAnnotation.from_data(y).summary)
                except TypeError:
                    print json.dumps(y)

            debug("roundtripped " + str(k) + " times")
            identical = Tests.compare_json(x,y)
            if not identical:
                debug("x:\n" + json.dumps(x))
                debug("y:\n" + json.dumps(y))

            # ensure that the result is identical
            self.failUnless(identical)
    
    @staticmethod
    def compare_json(x,y):
        try:
            assert type(x) == type(y)
            if isinstance(x,str) or isinstance(x,int) or isinstance(x,float) or isinstance(x,unicode):
                assert(x == y)
            elif isinstance(x,dict):
                assert x.keys() == y.keys()
                for i in x:
                    Tests.compare_json(x[i],y[i])
            else: # list
                assert len(x) == len(y)
                x.sort()
                y.sort()
                for i in range(len(x)):
                    Tests.compare_json(x[i],y[i])
        except AssertionError:
            return False
        return True
            
    
    def test_random_annotations_on_tree(self):

        ntax = random.randrange(2,100)
        
        for i in range(10):
            debug("generating annotations from tree " + str(i)) 
            tree = dendropy.simulate.birth_death_tree(birth_rate=0.1, death_rate=0, ntax=ntax)
            tree_label = str(i)

            # for each internal node, generate a random annotation
            annotations = []
            for j, n in enumerate(tree.internal_nodes()):
                n.label = j

                # generate the annotation, id matches the id applied to the node
                a = RandomAnnotation(id=j)
        
                # OVERWRITE the random target: include all ingroup, exclude all outgroup
                a.target = ReferenceTarget.from_data({
                    "type": "node",
                    "included_ids": [l.taxon.label for l in n.leaf_nodes()]
                })
                node_leafset = set(n.leaf_nodes())
                excluded = set()
                for l in tree.leaf_nodes():
                    if l not in node_leafset:
                        excluded.add(l)
                a.target.add_error_condition(TargetExcludesCondition(*[s.taxon.label for s in excluded]))

                annotations.append(a)

            # write the tree to a file
            test_tree_file = "tests/" + tree_label + ".input.tre"
            with open(test_tree_file, "w") as treefile:
                treefile.write(tree.as_string("newick"))
 
            # write the annotations to a file
            test_annotations_file = "tests/" + tree_label + ".input.annotations.json"
            with open(test_annotations_file, "w") as datafile:
                datafile.write("[\n")
                first = True
                for a in annotations:
                    if not first:
                        datafile.write(",\n")
                    else:
                        first = False
                    datafile.write(a.summary)
                datafile.write("\n]")

            # call the main method to process the data
            ot = open("tests/" + tree_label + ".output.tre", "w")
            ob = open("tests/" + tree_label + ".output.table", "w")
            main(test_tree_file, test_annotations_file, ot, ob, False)
     
            # check to see that the annotations were applied to the correct nodes (somehow)
            self.failUnless(True)
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('demo of muriqui mapping annotations to a tree')
    parser.add_argument('--taxon-tree',
                        type=int,
                        help='ott ID that will be used to fetch taxonomy/subtree to use as the tree to be annotated')
    parser.add_argument('--tree-node',
                        type=int,
                        help='treemachine node ID that will be used to fetch tree_of_life/subtree to use as the tree to be annotated')
    parser.add_argument('--tree-ott',
                        type=int,
                        help='ott ID that will be used to fetch tree_of_life/subtree to use as the tree to be annotated')
    parser.add_argument('--tree-file',
                        help='filepath to newick file with labels as ott IDs or using the name_ott#### convention')
    parser.add_argument('--out-table',
                        required=True,
                        help='file to output with the annotation placements')
    parser.add_argument('--out-tree',
                        required=True,
                        help='file to output with a tree with IDs to be used with the out-table')
    parser.add_argument('json', help='filepath to JSON file with annotations')
    args = parser.parse_args()
    annotations_file = args.json
    o_tree = open(args.out_tree, 'w')
    o_table = open(args.out_table, 'w')

    if args.tree_file is not None:
        tree_file = args.tree_file
    else:
        if args.taxon_tree is not None:
            _resp = TAXOMACHINE.subtree(int(args.taxon_tree))['subtree']
        elif args.tree_node is not None:
            _resp = TREEMACHINE.subtree(node_id=int(args.tree_node))['newick']
        elif args.tree_ott is not None:
            _resp = TREEMACHINE.subtree(ott_id=int(args.tree_ott))['newick']
        else:
            sys.exit('must specify a tree\n')
        import tempfile
        handle, tmpf = tempfile.mkstemp(suffix='.tre')
        sys.stderr.write('Writing fetched tree to "{}"'.format(tmpf))
        os.close(handle)
        _o = codecs.open(tmpf, 'w', encoding='utf-8')
        _o.write(_resp)
        _o.write(';\n')
        _o.close()
        tree_file = tmpf
    
    main(tree_file, annotations_file, o_tree, o_table)
