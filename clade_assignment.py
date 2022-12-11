import json, argparse
import urllib.request, json
from collections import defaultdict
import numpy as np
from config.weights import weights

def assign_alive(n, max_date, min_date):
    '''
    Assigns a node attribute "alive" to 1 if the node is within the date rate, 0 otherwise.
    '''
    if "children" in n:
        for c in n["children"]:
            assign_alive(c, max_date, min_date)
    else:
        numdate = n['node_attrs']['num_date']['value']
        n['alive'] = 1 if numdate<max_date and numdate>min_date else 0

def gather_tip_counts(n, attenuate=None, max_value=0, ignore_clade_breakpoints=False):
    '''
    Assigns a node attribute "ntips" to each node that counts the number of downstream tips.
    This can be attenuated similar to the LBI calculation.
    '''
    if "children" in n:
        tmp = 0
        for c in n["children"]:
            tmp_max = gather_tip_counts(c, attenuate=attenuate, max_value=max_value, ignore_clade_breakpoints=ignore_clade_breakpoints)
            if (ignore_clade_breakpoints and 'clade' in c['branch_attrs'].get('labels', {}))==False:
                tmp += c["ntips"]*attenuate(c)
            if tmp_max>max_value:
                max_value = tmp_max
        n["ntips"] = tmp
    else:
        n["ntips"] = 1 if n['alive'] else 0

    max_value = n["ntips"] if n["ntips"]>max_value else max_value

    return max_value

def assign_divergence(n, gene='HA1'):
    '''
    Assigns a node attribute "div" to each node that counts the number of mutations since the root of the tree.
    '''
    if "children" in n:
        for c in n["children"]:
            c['div'] = n['div'] + len(c['branch_attrs']['mutations'].get(gene,[]))
            assign_divergence(c)

def maxval(n, attr, val=-np.inf, best_node=None):
    if best_node is None:
        best_node=n
        val = n[attr]
    if "children" in n:
        for c in n["children"]:
            val, best_node = maxval(c, attr, val, best_node)

    better = n[attr]>val
    return n[attr] if better else val, n if better else best_node

def first_above(n, attr, cutoff):
    '''
    Returns the first node above the cutoff value for the given attribute when traversing the tree in pre-order.
    '''
    if n['node_attrs'][attr]['value']>cutoff:
        return n
    else:
        if "children" in n:
            for c in n["children"]:
                choice = first_above(c, attr, cutoff)
                if choice is not None:
                    return choice
        return None

def remove_attr(n, attr):
    if "children" in n:
        for c in n["children"]:
            remove_attr(c, attr)
    if attr in n:
        n.pop(attr)

def rename_clades(n, name_map):
    n['node_attrs']['trial_clade']['value'] = name_map(n['node_attrs']['trial_clade']['value'])
    if 'children' in n:
        for c in n['children']:
            rename_clades(c, name_map)

def rename_clade_labels(n, name_map):
    if "labels" in n['branch_attrs'] and 'trial_clade' in n['branch_attrs']['labels']:
        print(n['branch_attrs']['labels']['trial_clade'])
        n['branch_attrs']['labels']['trial_clade'] = name_map(n['branch_attrs']['labels']['trial_clade'])
    if 'children' in n:
        for c in n['children']:
            rename_clade_labels(c, name_map)

def score(n, weights=None, ntip_scale=1):
    if weights is None:
        weights = {}
    score = np.sqrt(n['ntips']/ntip_scale)
    mut_weight = sum([weights.get(int(x[1:-1]),1) for x in n['branch_attrs']['mutations'].get('HA1',[])])
    if mut_weight>0:
        score += mut_weight/(4 + mut_weight)
    else:
        score = 0.0
    return score


def assign_score(n, score=None, weights=None, ntip_scale=1):
    if "children" in n:
        for c in n["children"]:
            assign_score(c, score=score, weights=weights, ntip_scale=ntip_scale)
    n['node_attrs']["score"] = {'value': score(n, weights=weights, ntip_scale=ntip_scale)}


def assign_clade(n, clade):
    if "children" in n:
        for c in n["children"]:
            assign_clade(c, clade)
    n['node_attrs']["trial_clade"] = {'value': clade}

def assign_clades_to_branch(n, hierarchy, cutoff=1.0):
    node = first_above(n, 'score', cutoff)
    if node is not None:
        if "clade_score" not in node["node_attrs"]:
            print("new clade",node['name'], node['div'], node['ntips'], node['node_attrs']['score']['value'])

            if 'labels' not in node['branch_attrs']:
                node['branch_attrs']['labels'] = {}

            node["node_attrs"]["clade_score"] = {"value":node["node_attrs"]["score"]["value"]}
            print(node["node_attrs"]["clade_score"])

            parent_clade = node['node_attrs']["trial_clade"]["value"]
            new_suffix = len(hierarchy[parent_clade])+1
            new_clade = tuple(list(parent_clade) + [new_suffix])
            hierarchy[parent_clade].append(new_suffix)
            node['branch_attrs']['labels']['trial_clade'] = new_clade
            assign_clade(node, new_clade)

            hierarchy[new_clade] = []
            clades[new_clade] = node

        if 'children' in n:
            for c in n["children"]:
                hierarchy = assign_clades_to_branch(c, hierarchy, cutoff=cutoff)

    return hierarchy

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Assign clades to a tree")
    parser.add_argument('--lineage', default='h3n2', help="lineage to assign clades to")
    parser.add_argument('--resolution', default='6y')
    parser.add_argument('--output', default='annotated_tree.json')

    args = parser.parse_args()

    with urllib.request.urlopen(f"https://nextstrain.org/charon/getDataset?prefix=flu/seasonal/{args.lineage}/ha/{args.resolution}") as url:
        data = json.loads(url.read())

    T = data["tree"]

    max_date, min_date = 2030,2010

    clades = {}
    assign_alive(T,max_date,min_date)
    max_value = gather_tip_counts(T, lambda x:np.exp(-len(x['branch_attrs']['mutations'].get('HA1',[]))/1), ignore_clade_breakpoints=True)

    root_clade = ('A',)
    T['div']=0
    assign_divergence(T)
    assign_score(T, score, weights=weights[args.lineage], ntip_scale=max_value)




    assign_clade(T, root_clade)
    hierarchy = {root_clade:[]}
    hierarchy = assign_clades_to_branch(T, hierarchy, cutoff=0.75)

    aliases = {root_clade: "R"}
    major_clades = defaultdict(list)
    alpha = 'ABCDEFGHIJKLMNOPQRST'
    for clade in clades:
        node  = clades[clade]
        print(node["node_attrs"]["clade_score"]["value"])
        if node["node_attrs"]["clade_score"]["value"]>1:
            year = int(node["node_attrs"]["num_date"]["value"]) - 2000
            letter = alpha[len(major_clades[year])%len(alpha)]
            aliases[node["node_attrs"]["trial_clade"]["value"]] = f"{year}{letter}"
            major_clades[year].append(letter)


    #clade_map = list(range(clade_counter))
    #random.shuffle(clade_map)
    name_map = {}
    for clade in hierarchy:
        if clade in aliases:
            name_map[clade] = aliases[clade]
        else:
            prefix_length = -1
            for alias in aliases:
                if len(alias)>prefix_length:
                    if clade[:len(alias)] == alias:
                        best_alias = alias
                        prefix_length = len(alias)
            name_map[clade] = aliases[best_alias] + '.' + '.'.join([str(x) for x in clade[prefix_length:]])

    rename_clades(T, lambda x:name_map[x])
    rename_clade_labels(T, lambda x:name_map[x])

    data['meta']['colorings'].append({'key':"trial_clade", 'type':'ordinal', 'title':"trial clade"})
    data['meta']['colorings'].append({'key':"score", 'type':'continuous', 'title':"clade score"})
    with open(f'auspice/trial_{args.lineage}_clades.json', 'w') as fh:
        json.dump(data, fh)
