import argparse
import random
import sys
import dendropy

def refine(T, t):
    T2 = T.clone()    
    T2.retain_taxa([i.taxon for i in t.leaf_nodes()])

    splits = set([i for i in t.encode_bipartitions()])

    

    for b in T2.encode_bipartitions():
        if t.is_compatible_with_bipartition(b) and b not in splits:
            splits.add(b)

    
    splits = [i.split_bitmask for i in splits]

    tout = dendropy.Tree.from_split_bitmasks(splits, t.taxon_namespace)
    tout.retain_taxa([i.taxon for i in t.leaf_nodes()])
    return tout
            
def edges_by_split(t, split):
    
    taxonset = t.seed_node.leafset_bitmask

    return [e for e in t.edges() 
            if (e.split_bitmask & taxonset) == split 
            or ((taxonset ^ e.split_bitmask) & taxonset) == split]

def reconcile(T, t, random_path):
    assert(len(T.leaf_nodes()) == len(t.leaf_nodes()) + 1)
    T.encode_bipartitions()
    t.encode_bipartitions()

    tleafs = set([i.taxon for i in t.leaf_nodes()])
    n = [i for i in T.leaf_nodes() if i.taxon not in tleafs][0]

    T.reroot_at_edge(n.edge)
    T.reroot_at_edge(n.edge)

    TA = dendropy.Tree(T)
    TA.retain_taxa([i.taxon for i in t.leaf_nodes()])
    TA.encode_bipartitions()

    t_taxa = set([i.taxon for i in t.leaf_nodes()])
 
    if random_path:
        root_node = TA.seed_node.child_nodes()[random.randint(0, 1)]
    else:
        root_node = TA.seed_node.child_nodes()[0]
    root_edge = root_node.edge
    root_split = root_edge.split_bitmask
    root_split_t = edges_by_split(t, root_split)
    while(len(root_split_t) == 0):
        if random_path:
            root_node = root_node.child_nodes()[random.randint(0, 1)]
        else:
            root_node = root_node.child_nodes()[0]
        root_edge = root_node.edge
        root_split = root_edge.split_bitmask
        root_split_t = edges_by_split(t, root_split)

    t.reroot_at_edge(root_split_t[0])
    t.reroot_at_edge(root_split_t[0])

    t.seed_node.add_child(n)
    return t

def complete(T, t, random_path):
    tleafs = set([i.taxon for i in t.leaf_nodes()])
    Tleafs = set([i.taxon for i in T.leaf_nodes()])

    if tleafs == Tleafs:
        return t
    
    n = [i for i in T.leaf_nodes() if i.taxon not in tleafs][0]

    Tx = dendropy.Tree(T)
    Tx.prune_taxa([n.taxon])
    Tx.update_bipartitions()

    t = complete(Tx, t, random_path)
    t = reconcile(T, t, random_path)

    return t


def main(args):
    output_file = open(args.output, 'w')
    tns = dendropy.TaxonNamespace()
    trees = dendropy.TreeList.get(path=args.input, schema='newick', taxon_namespace=tns)
    T = dendropy.Tree.get(path=args.backbone, schema='newick', taxon_namespace=tns)
    
    for t in trees:
        if args.refine:
            t = refine(dendropy.Tree(T), dendropy.Tree(t))
        t.resolve_polytomies(limit=2, update_bipartitions=True)
        tx = complete(dendropy.Tree(T), dendropy.Tree(t), args.random_path)
        output_file.write(str(tx) + ';' + '\n')

    output_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
                        help='File containing incomplete trees',
                        required=True)
    parser.add_argument('-o','--output', type=str,
                        help='Name of output file',
                        required=True)
    parser.add_argument('-b','--backbone', type=str,
                        help='File containing the backbone tree',
                        required=True)
    parser.add_argument('-r','--random_path',
                        action="store_true", default=False,
                        help='Complete tree using a random path')
    parser.add_argument('-f','--refine',
                        action="store_true", default=False,
                        help='Refine tree before completing')
    args = parser.parse_args()
    main(args)
