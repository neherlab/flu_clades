import json, argparse
from collections import defaultdict
import numpy as np


def prepare_tree(n, max_date, min_date):
    '''
    Assigns a node attribute "alive" to 1 if the node is within the date rate, 0 otherwise.
    '''
    if "children" in n:
        tmp_alive = False
        for c in n["children"]:
            prepare_tree(c, max_date, min_date)
            tmp_alive = tmp_alive or c['alive']
        n['alive'] = tmp_alive
    else:
        try:
            numdate = n['node_attrs']['num_date']['value']
            n['alive'] = 1 if numdate<max_date and numdate>min_date else 0
        except:
            n['alive'] = 1

    n['clade_break_point'] = False


def label_backbone(tree, key):
    '''
    This function labels all branches/nodes from the root to existing clade labels
    as 'back-bone'. This can be used to prevent introduction of additional clades
    along the backbone
    '''
    def label_backbone_recursive(n, key):
        on_backbone = False
        if "children" in n:
            for c in n["children"]:
                label_backbone_recursive(c, key)
                if c["backbone"]:
                    on_backbone = True

        if on_backbone or n["clade_break_point"]:
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


def assign_divergence(n, genes=['HA1']):
    '''
    Assigns a node attribute "div" to each node that counts the number of mutations
    since the root of the tree
    '''
    if "children" in n:
        for c in n["children"]:
            c['div'] = n['div']
            for gene in genes:
                bad_states = ['-', 'N'] if gene=='nuc' else ['X', '-']
                c['div'] += len([x for x in c['branch_attrs']['mutations'].get(gene,[])
                                if x[0] not in bad_states or x[-1] not in bad_states])
            assign_divergence(c, genes)


def calc_phylo_score(n, distance=None, ignore_backbone=False):
    '''
    Assigns a node attribute "bushiness" to each node that counts the number of downstream tips.
    This can be distanced similar to the LBI calculation.
    '''
    if "children" in n:
        n["bushiness"] = 0
        n["ntips"] = 0
        for c in n["children"]:
            calc_phylo_score(c, distance=distance, ignore_backbone=ignore_backbone)

            if (ignore_backbone and c["backbone"])==False:
                n["bushiness"] += c["bushiness"]*np.exp(-distance(c))
                n["bushiness"] += (1-np.exp(-distance(c))) if c['alive'] else 0

            n["ntips"] += c["ntips"]
    else:
        n["bushiness"] = 1 if n['alive'] else 0
        n["ntips"] = 1

    n['node_attrs']["bushiness_raw"] = {'value': n["bushiness"]}


def calc_phylo_scale(T):
    '''
    Calculates the phylogenetic scale of the tree as the median bushiness of all internal nodes
    '''
    def collect_recursive(n, values):
        if "children" in n and len(n["children"])>0:
            for c in n["children"]:
                collect_recursive(c, values)
            if n['alive']: # only include alive nodes and skip terminals
                values.append(n["bushiness"])

    values = []
    collect_recursive(T, values)
    from scipy.stats import scoreatpercentile
    return scoreatpercentile(values,80)


def score(n, weights=None, bushiness_scale=1, ignore_backbone=False,
          proteins=None, branch_length_scale=4):
    '''
    Assign a score to each node that combines the phylogenetic score and the branch length
    '''
    if weights is None:
        weights = {}
    if proteins is None:
        proteins = []

    score = n["bushiness"]/(n["bushiness"] + bushiness_scale)
    n['node_attrs']['bushiness'] = {'value': score}

    mut_weight = 0
    for cds in proteins:
        w = weights.get(cds, {})
        for x in n['branch_attrs']['mutations'].get(cds,[]):
            if x[-1] in ['-', 'X']: continue
            pos = int(x[1:-1])
            mut_weight += w[pos] if pos in w else w.get("default", 0)

    n['node_attrs']['branch_score'] = {'value': mut_weight/(branch_length_scale + mut_weight)}
    score += mut_weight/(branch_length_scale + mut_weight)

    # return 0 if the node is on the backbone and we are ignoring backbone nodes
    if ignore_backbone and n['backbone']:
        return 0.0
    # return 0 if there are no mutations in the proteins of interest
    if sum([len(n['branch_attrs']['mutations'].get(cds,[])) for cds in proteins])==0:
        return 0.0

    return score


def assign_score(n, score=None, **kwargs):
    '''
    recursively assign the clade demarcation score to each branch
    '''
    if "children" in n:
        for c in n["children"]:
            assign_score(c, score=score, **kwargs)
    n['node_attrs']["score"] = {'value': score(n, **kwargs)}


def assign_clade(n, clade, key):
    '''Assign a clade to a node and recursively to all its children'''
    if "children" in n:
        for c in n["children"]:
            assign_clade(c, clade, key)
    n['node_attrs'][key] = {'value': clade}


def assign_new_clades_to_branches(n, hierarchy, new_key, new_clades=None,
                                  cutoff=1.0, divergence_addition=None, divergence_base=0,
                                  divergence_scale=4, min_size=5):
    '''
    walk through the tree in pre-order (by recursively calling this function)
    and call a new clade whenever there is a branch that crosses the threshold
    '''
    if divergence_addition:
        delta_div = n['div']-divergence_base  # calculate the divergence since the parent
        div_score = divergence_addition*delta_div/(delta_div+divergence_scale)
    else: div_score=0
    n["node_attrs"]['div_score'] = {'value': div_score}

    if (n["node_attrs"]['score']['value'] + div_score > cutoff) and (n["ntips"]>min_size):
        if 'labels' not in n['branch_attrs']:
            n['branch_attrs']['labels'] = {}

        # determine parent clade
        parent_clade = tuple(n['node_attrs'][f"full_{new_key}"]["value"])
        if parent_clade in hierarchy:
            # determine the number of existing children of the parent and the index of the new subclade
            new_suffix = len(hierarchy[parent_clade])+1
            if new_suffix>2: # consistency check
                assert new_suffix==hierarchy[parent_clade][-1]+1

            hierarchy[parent_clade].append(new_suffix)
            new_clade = tuple(list(parent_clade) + [new_suffix])
        else:
            new_clade = parent_clade
            print("new clade found, but no parent to attach it to")

        hierarchy[new_clade] = []
        new_clades[new_clade] = n
        assign_clade(n, new_clade, f"full_{new_key}")
        n["clade_break_point"] = True  # mark as clade break_point

    # reset divergence to clade break point.
    if n['clade_break_point']:
        # print(n['div'] - divergence_base, n['branch_attrs'].get('labels',{}))
        divergence_base=n['div']

    if 'children' in n:
        for c in n["children"]:
            assign_new_clades_to_branches(c, hierarchy, new_key,
                                new_clades=new_clades, cutoff=cutoff,
                                divergence_addition=divergence_addition,
                                divergence_base=divergence_base,
                                divergence_scale=divergence_scale,
                                min_size=min_size)


def copy_over_old_clades(tree, old_key, new_key):
    def copy_recursive(n,old_key, new_key):
        n["node_attrs"][new_key] = {k:v for k,v in n["node_attrs"][old_key].items()}
        n["node_attrs"][f"full_{new_key}"] = {k:v for k,v in n["node_attrs"][f"full_{old_key}"].items()}
        n['clade_break_point'] = False
        if old_key in n["branch_attrs"].get("labels",{}):
            n["branch_attrs"]["labels"][new_key] = n["branch_attrs"]["labels"][old_key]
            n['clade_break_point'] = True

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

def get_clade_map(fname):
    if fname is None:
        return {}, {}

    aliases = {}
    with open(fname) as fh:
        old_to_new_clades = json.load(fh)
    for k,v in old_to_new_clades.items():
        if len(v[0])==1:
            aliases[tuple(v[1])] = v[0]
    short_to_full_clades = {v[0]:v[1] for v in old_to_new_clades.values()}

    return short_to_full_clades, aliases

def get_tree(fname, max_date=None, min_date=None, add_to_existing=False, old_key=None, new_key=None):
    with open(fname) as fh:
        data = json.load(fh)

    T = data["tree"]
    prepare_tree(T, max_date=max_date, min_date=min_date)
    T['div']=0
    assign_divergence(T, ['nuc'])

    hierarchy = defaultdict(list)
    if add_to_existing:
        print("adding existing clades")
        existing_clades, existing_full_clades = get_existing_clade_labels(T, old_key)
        for full_clade in existing_full_clades:
            if len(full_clade)>1:
                hierarchy[full_clade[:-1]].append(full_clade[-1])
            if full_clade not in hierarchy:
                hierarchy[full_clade] = []
        copy_over_old_clades(T, old_key, new_key)
        label_backbone(T, old_key)
    else:
        hierarchy[('A',)] = ['A']
        print("skipping existing clades")
        # set root to be clade A
        T['node_attrs'][f'full_{args.new_key}'] = hierarchy[('A',)]
        assign_clade(T, hierarchy[('A',)], f"full_{new_key}")
        assign_clade(T, 'A', new_key)

    hierarchy = {k:sorted(hierarchy[k]) for k in sorted(hierarchy.keys())}

    return data, T, hierarchy

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Assign clades to a tree")
    parser.add_argument('--input', help="json to assign clades to")
    parser.add_argument('--lineage', default='h3n2', help="json to assign clades to")
    parser.add_argument('--segment', default='ha', help="json to assign clades to")
    parser.add_argument('--add-to-existing', action='store_true', help="respect old clades")
    parser.add_argument('--clade-map',help="json with renamed clades")
    parser.add_argument('--weights',help="json with mutation weights")
    parser.add_argument("--old-key", type=str, help='name to save clades under')
    parser.add_argument("--new-key", type=str, help='name to save clades under')
    parser.add_argument('--output', default='recladed_tree.json')

    args = parser.parse_args()

    if 'rsv' in args.lineage.lower():
        proteins = ['G', 'F', 'L', 'N', 'P', 'M']
        max_date, min_date = 2030,1990
        cutoff=1.3
        bushiness_branch_scale = 10.0
        divergence_scale = 20
        branch_length_scale = 8
        divergence_addition = 1.0
        min_size = 10
    else:
        proteins = ['HA1', 'HA2'] if args.segment == 'ha' else ['NA']
        max_date, min_date = 2030,2020
        cutoff=1.0
        bushiness_branch_scale = 1
        divergence_scale = 8
        branch_length_scale = 4
        divergence_addition = 1.0
        min_size = 20

    if args.add_to_existing:
        short_to_full_clades, aliases = get_clade_map(args.clade_map)
    else:
        short_to_full_clades, aliases = {}, {}

    with open(args.weights) as fh:
        weights = json.load(fh)

    data, T, hierarchy = get_tree(args.input, max_date=max_date, min_date=min_date,
                            add_to_existing=args.add_to_existing,
                            old_key=args.old_key, new_key=args.new_key)


    # nucleotide branch length excluding gaps and N
    branch_length_function = lambda x:len([y for y in x['branch_attrs']['mutations'].get('nuc',[])
                                           if y[-1] not in ['N', '-'] and y[0] not in ['N', '-']])/bushiness_branch_scale
    calc_phylo_score(T, branch_length_function, ignore_backbone=args.add_to_existing)
    bushiness_scale = calc_phylo_scale(T)
    print("phylo_score_scale", bushiness_scale)

    # compute aggregate score from branches and phylo/bushiness
    assign_score(T, score, weights=weights[args.lineage],
                 bushiness_scale=bushiness_scale, ignore_backbone=args.add_to_existing,
                 proteins=proteins, branch_length_scale=branch_length_scale)

    # assign clades while also taking into account the divergence, modifies hierarchy and new_clades in place
    new_clades = {}
    assign_new_clades_to_branches(T, hierarchy, args.new_key,
        new_clades=new_clades, cutoff=cutoff, divergence_addition=divergence_addition,
        divergence_base=0.0, divergence_scale=divergence_scale, min_size=min_size)

    # process and assign human readable clade names
    for new_clade in new_clades:
        clade_name = full_clade_to_short_name(new_clade, aliases)
        n = new_clades[new_clade]
        if "labels" not in n["branch_attrs"]: n["branch_attrs"]["labels"] = {}
        n["branch_attrs"]["labels"][args.new_key] = clade_name
        assign_clade(n, clade_name, args.new_key)
        print("suggested clade:", clade_name,
              {k:v for k, v in n["branch_attrs"]["mutations"].items() if k!='nuc'})

    # export
    data['meta']['colorings'].append({'key':args.new_key, 'type':'ordinal', 'title':args.new_key})
    data['meta']['colorings'].append({'key':"score", 'type':'continuous', 'title':"clade score"})
    data['meta']['colorings'].append({'key':"bushiness_raw", 'type':'continuous', 'title':"phylo_score_raw"})
    data['meta']['colorings'].append({'key':"bushiness", 'type':'continuous', 'title':"phylo_score"})
    data['meta']['colorings'].append({'key':"div_score", 'type':'continuous', 'title':"div_score"})
    data['meta']['colorings'].append({'key':"branch_score", 'type':'continuous', 'title':"branch_score"})
    data['meta']['filters'].append('new_clade')

    with open(args.output, 'w') as fh:
        json.dump(data, fh, indent=0)
