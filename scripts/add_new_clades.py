import json, argparse
from collections import defaultdict
import numpy as np


def label_backbone(tree, key):
    def label_backbone_recursive(n, key):
        on_backbone = False
        if "children" in n:
            for c in n["children"]:
                label_backbone_recursive(c, key)
                if c["backbone"]:
                    on_backbone = True

        if on_backbone or ("labels" in n["branch_attrs"] and key in n["branch_attrs"]["labels"]):
            n["backbone"] = True
        else:
            n["backbone"] = False

    label_backbone_recursive(tree, key)

def get_existing_clade_labels(tree, key):
    '''
    Returns a list of existing clade labels in the tree.
    '''
    def add_clades(n, clades, full_clades, key):
        if "labels" in n["branch_attrs"] and key in n["branch_attrs"]["labels"]:
            clades.add(n["branch_attrs"]["labels"][key])
        if f"full_{key}" in n["node_attrs"]:
            full_clades.add(tuple(n["node_attrs"][f"full_{key}"]["value"]))

        if "children" in n:
            for c in n["children"]:
                add_clades(c, clades, full_clades, key)

    clades = set()
    full_clades = set()
    add_clades(tree, clades, full_clades, key)
    return clades, full_clades

def assign_divergence(n, genes=['HA1'], key=None):
    '''
    Assigns a node attribute "div" to each node that counts the number of mutations since the root of the tree.
    '''
    if key and key in n["branch_attrs"].get("labels",{}):
        n['div'] = 0
    if "children" in n:
        for c in n["children"]:
            c['div'] = n['div'] + sum([
                  len([x for x in c['branch_attrs']['mutations'].get(gene,[]) if x[-1] not in ['-', 'X']])
                for gene in genes])
            assign_divergence(c, genes, key)


def assign_alive(n, max_date, min_date):
    '''
    Assigns a node attribute "alive" to 1 if the node is within the date rate, 0 otherwise.
    '''
    if "children" in n:
        for c in n["children"]:
            assign_alive(c, max_date, min_date)
    else:
        try:
            numdate = n['node_attrs']['num_date']['value']
        except:
            numdate = 2022
        n['alive'] = 1 if numdate<max_date and numdate>min_date else 0


def gather_tip_counts(n, distance=None, scale=1, max_value=0, ignore_backbone=False):
    '''
    Assigns a node attribute "bushiness" to each node that counts the number of downstream tips.
    This can be distanced similar to the LBI calculation.
    '''
    if "children" in n:
        n["bushiness"] = 0
        n["ntips"] = 0
        for c in n["children"]:
            tmp_max = gather_tip_counts(c, distance=distance, scale=scale,
                        max_value=max_value, ignore_backbone=ignore_backbone)
            if tmp_max>max_value:
                max_value = tmp_max

            if (ignore_backbone and c["backbone"])==False:
                n["bushiness"] += c["bushiness"]*np.exp(-distance(c)/scale)
                n["bushiness"] += (1-np.exp(-distance(c)/scale))
            n["ntips"] += c["ntips"]
    else:
        n["bushiness"] = (1-np.exp(-distance(n)/scale)) if n['alive'] else 0
        n["ntips"] = 1

    n['node_attrs']["phylo"] = {'value': n["bushiness"]}
    max_value = n["bushiness"] if n["bushiness"]>max_value else max_value

    return max_value

def score(n, weights=None, bushiness_scale=1, ignore_backbone=False,
          genes=None, core_genes=None, branch_length_scale=4):
    if weights is None:
        weights = {}
    if genes is None:
        genes = []
    if core_genes is None:
        core_genes = []

    n['tip_score'] = np.sqrt(n["bushiness"]/bushiness_scale)
    if ignore_backbone and n['backbone']:
        return 0.0
    elif sum([len(n['branch_attrs']['mutations'].get(gene,[])) for gene in genes])==0:
        return 0.0

    score = n['tip_score']
    mut_weight = 0
    for gene in core_genes:
        mut_weight += sum([weights.get(int(x[1:-1]),1)
                    for x in n['branch_attrs']['mutations'].get(gene,[]) if x[-1] not in ['-', 'X']
                    ])

    score += mut_weight/(branch_length_scale + mut_weight)
    return score


def assign_score(n, score=None, **kwargs):
    if "children" in n:
        for c in n["children"]:
            assign_score(c, score=score, **kwargs)
    n['node_attrs']["score"] = {'value': score(n, **kwargs)}

def assign_clade(n, clade, key):
    if "children" in n:
        for c in n["children"]:
            assign_clade(c, clade, key)
    n['node_attrs'][key] = {'value': clade}

def assign_new_clades_to_branches(n, hierarchy, new_key, new_clades=None,
                                  cutoff=1.0, divergence_addition=None, divergence_base=0,
                                  divergence_scale=4, min_size=5):
    if divergence_addition:
        div_score = divergence_addition*(n['div']-divergence_base)/((n['div']-divergence_base)+divergence_scale)*n['tip_score']
    else: div_score=0
    #print(n["node_attrs"]['score']['value'],  div_score, cutoff)
    if n["node_attrs"]['score']['value'] + div_score > cutoff and n["ntips"]>min_size:
        print('thres')
        divergence_base=n['div']
        if 'labels' not in n['branch_attrs']:
            n['branch_attrs']['labels'] = {}

        parent_clade = tuple(n['node_attrs'][f"full_{new_key}"]["value"])
        print(parent_clade)
        if parent_clade in hierarchy:
            new_suffix = len(hierarchy[parent_clade])+1
            print(hierarchy[parent_clade], new_suffix)
            if new_suffix>2:
                assert new_suffix==hierarchy[parent_clade][-1]+1

            new_clade = tuple(list(parent_clade) + [new_suffix])
            hierarchy[parent_clade].append(new_suffix)
            hierarchy[new_clade] = []
            new_clades[new_clade] = n

            assign_clade(n, new_clade, f"full_{new_key}")

    if 'children' in n:
        for c in n["children"]:
            hierarchy = assign_new_clades_to_branches(c, hierarchy, new_key,
                                new_clades=new_clades, cutoff=cutoff,
                                divergence_addition=divergence_addition,
                                divergence_base=divergence_base)

    return hierarchy

def copy_over_old_clades(tree, old_key, new_key):
    def copy_recursive(n,old_key, new_key):
        n["node_attrs"][new_key] = {k:v for k,v in n["node_attrs"][old_key].items()}
        n["node_attrs"][f"full_{new_key}"] = {k:v for k,v in n["node_attrs"][f"full_{old_key}"].items()}
        if old_key in n["branch_attrs"].get("labels",{}):
            n["branch_attrs"]["labels"][new_key] = n["branch_attrs"]["labels"][old_key]

        if "children" in n:
            for c in n["children"]:
                copy_recursive(c, old_key, new_key)

    copy_recursive(tree, old_key, new_key)

def full_clade_to_short_name(full_clade, aliases):
    for full_parent in sorted(aliases.keys(), key=lambda x:len(x), reverse=True):
        if tuple(full_clade[:len(full_parent)])==full_parent:
            clade_name = aliases[full_parent]
            if len(full_parent)<len(full_clade):
                clade_name += '.' + '.'.join([str(x) for x in full_clade[len(full_parent):]])
            return clade_name
    return '.'.join([str(x) for x in full_clade])

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Assign clades to a tree")
    parser.add_argument('--input', help="json to assign clades to")
    parser.add_argument('--lineage', default='h3n2', help="json to assign clades to")
    parser.add_argument('--add-to-existing', action='store_true', help="respect old clades")
    parser.add_argument('--clade-map',help="json with renamed clades")
    parser.add_argument('--weights',help="json with mutation weights")
    parser.add_argument("--old-key", type=str, help='name to save clades under')
    parser.add_argument("--new-key", type=str, help='name to save clades under')
    parser.add_argument('--output', default='recladed_tree.json')

    args = parser.parse_args()

    if 'rsv' in args.lineage.lower():
        all_genes = ['G', 'F', 'L', 'N', 'P', 'M']
        core_genes = ['G', 'F']
        tip_count_divergence_scale = 10.0
        divergence_addition=0.7
        max_date, min_date = 2030,1990
        cutoff=1.0
        divergence_scale=20
        branch_length_scale = 4
        min_size = 10
    else:
        all_genes = ['HA1', 'HA2']
        core_genes = ['HA1']
        tip_count_divergence_scale = 1
        max_date, min_date = 2030,2020
        cutoff=0.85
        divergence_addition=1.0
        divergence_scale=4
        branch_length_scale = 4
        min_size = 5


    with open(args.clade_map) as fh:
        old_to_new_clades = json.load(fh)
    with open(args.weights) as fh:
        weights = json.load(fh)

    with open(args.input) as fh:
        data = json.load(fh)

    aliases = {}
    for k,v in old_to_new_clades.items():
        if len(v[0])==1:
            aliases[tuple(v[1])] = v[0]
    short_to_full_clades = {v[0]:v[1] for v in old_to_new_clades.values()}

    T = data["tree"]
    assign_alive(T, max_date=max_date, min_date=min_date)

    hierarchy = defaultdict(list)
    if args.add_to_existing:
        print("adding existing clades")
        existing_clades, existing_full_clades = get_existing_clade_labels(T, args.old_key)
        label_backbone(T, args.old_key)

        major_clades = sorted(set([x[0] for x in existing_clades]))
        for full_clade in existing_full_clades:
            if len(full_clade)<=1:
                continue
            hierarchy[full_clade[:-1]].append(full_clade[-1])
            if full_clade not in hierarchy:
                hierarchy[full_clade] = []
        copy_over_old_clades(T, args.old_key, args.new_key)
    else:
        hierarchy[('A',)] = ['A']
        print("skipping existing clades")
        T['node_attrs'][f'full_{args.new_key}'] = hierarchy[('A',)]
        assign_clade(T, hierarchy[('A',)], f"full_{args.new_key}")


    hierarchy = {k:sorted(hierarchy[k]) for k in sorted(hierarchy.keys())}
    T['div']=0
    assign_divergence(T, core_genes, args.old_key if args.add_to_existing else None)
    max_value = gather_tip_counts(T,
                    lambda x:len([y for y in x['branch_attrs']['mutations'].get('nuc',[]) if y[-1] not in ['N', '-'] and y[0] not in ['N', '-']]),
                    scale=tip_count_divergence_scale,
                    ignore_backbone=args.add_to_existing)
    print("phylo_score_max", max_value)
    assign_score(T, score, weights=weights[args.lineage],
                 bushiness_scale=max_value, ignore_backbone=args.add_to_existing,
                 genes=all_genes, core_genes=core_genes, branch_length_scale=branch_length_scale)

    new_clades = {}
    hierarchy = assign_new_clades_to_branches(T, hierarchy, args.new_key,
                    new_clades=new_clades, cutoff=cutoff, divergence_addition=divergence_addition,
                    divergence_base=0.0, divergence_scale=divergence_scale, min_size=min_size)

    for new_clade in new_clades:
        clade_name = full_clade_to_short_name(new_clade, aliases)
        n = new_clades[new_clade]
        if "labels" not in n["branch_attrs"]: n["branch_attrs"]["labels"] = {}
        n["branch_attrs"]["labels"][args.new_key] = clade_name
        assign_clade(n, clade_name, args.new_key)
        print("suggested clade:", clade_name, {k:v for k, v in n["branch_attrs"]["mutations"].items() if k!='nuc'})

    data['meta']['colorings'].append({'key':args.new_key, 'type':'ordinal', 'title':args.new_key})
    data['meta']['colorings'].append({'key':"score", 'type':'continuous', 'title':"clade score"})
    data['meta']['colorings'].append({'key':"phylo", 'type':'continuous', 'title':"phylo_weight"})
    data['meta']['filters'].append('new_clade')

    with open(args.output, 'w') as fh:
        json.dump(data, fh, indent=0)
