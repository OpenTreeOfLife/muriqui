import dendropy, random, argparse

import muriqui as m

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="generate random test trees and annotations")

    parser.add_argument("-l", "--label", nargs=1, required=True, help="The label to use for the output files.")

    args = parser.parse_args()
    
    label = args.label[0]

    # generate a random tree
    tree = dendropy.simulate.birth_death_tree(birth_rate=0.1,death_rate=0,ntax=8)

    # for each internal node, generate a random annotation
    annotations = []
    for i, n in enumerate(tree.internal_nodes()):
        n.label = i

        # generate the annotation
        a = m.PhyloReferencedAnnotation()
        a.id = n.label
        a.annotated_at = "datetime"
        a.annotated_by = m.Entity.from_data({"name": "cody"})
        a.body = {"test_value": random.random()}
        
        # set the target: include all ingroup, exclude all outgroup
        a.target = m.ReferenceTarget.from_data({
            "type": "node",
            "included_ids": [l.taxon.label for l in n.leaf_nodes()]
        })
        node_leafset = set(n.leaf_nodes())
        for l in tree.leaf_nodes():
            if l not in node_leafset:
                a.target.add_error_condition(m.CladeExcludesCheck(l.taxon.label))

        annotations.append(a)
    
    # write the tree to a file
    with open(label+".labeled.tre", "w") as treefile:
        treefile.write(tree.as_string("newick"))
    
    # write the annotations to a file
    with open(label+".annotations.json", "w") as datafile:        
        for a in annotations:
            datafile.write(a.summary + "\n")