import json, argparse
import urllib.request, json
from collections import defaultdict
import numpy as np


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

def gather_tip_counts(n, attenuate=None):
    '''
    Assigns a node attribute "ntips" to each node that counts the number of downstream tips.
    This can be attenuated similar to the LBI calculation.
    '''
    if "children" in n:
        tmp = 0
        for c in n["children"]:
            gather_tip_counts(c, attenuate=attenuate)
            tmp += c["ntips"]*attenuate(c)
        n["ntips"] = tmp
    else:
        n["ntips"] = 1 if n['alive'] else 0

def assign_divergence(n, gene='HA1'):
    '''
    Assigns a node attribute "div" to each node that counts the number of mutations since the root of the tree.
    '''
    if "children" in n:
        for c in n["children"]:
            c['div'] = n['div'] + len(c['branch_attrs']['mutations'].get('HA1',[]))
            assign_divergence(c)

def assign_score(n, score=None):
    if "children" in n:
        for c in n["children"]:
            assign_score(c, score=score)
    n['node_attrs']["score"] = {'value': score(n)}


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

def assign_clade(n, clade):
    if "children" in n:
        for c in n["children"]:
            assign_clade(c, clade)
    n['node_attrs']["trial_clade"] = {'value': clade}

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

def score(n):
    if n['div']<2:
        return 0.0

    score = n['ntips']
    nmuts = len(n['branch_attrs']['mutations'].get('HA1',[]))
    score += 20*nmuts/(4 + nmuts)
    # if 'labels' in n['branch_attrs']:
    #     if "clade" in n['branch_attrs']['labels']:
    #         score += 100
    return score

def assign_clades_to_branch(n, hierarchy):
    node = first_above(n, 'score', 20)
    if node is not None:
        print("new clade",node['name'], node['div'], node['ntips'], node['node_attrs']['score']['value'])
        node['div']=0
        if 'labels' not in node['branch_attrs']:
            node['branch_attrs']['labels'] = {}
        assign_divergence(node)
        node["node_attrs"]["clade_score"] = {"value":node["node_attrs"]["score"]["value"]}
        print(node["node_attrs"]["clade_score"])
        assign_score(node, score)
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
            hierarchy = assign_clades_to_branch(c, hierarchy)

        print(node["node_attrs"]["clade_score"])

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
    gather_tip_counts(T, lambda x:np.exp(-len(x['branch_attrs']['mutations'].get('HA1',[]))/1))

    root_clade = ('A',)
    T["node_attrs"]["clade_score"] = {"value": 100}
    T['div']=0
    assign_divergence(T)
    assign_score(T, score)
    assign_clade(T, root_clade)
    hierarchy = {root_clade:[]}
    hierarchy = assign_clades_to_branch(T, hierarchy)

    aliases = {root_clade: "R"}
    major_clades = defaultdict(list)
    alpha = 'ABCDEFGHIJKLMNOPQRST'
    for clade in clades:
        node  = clades[clade]
        print(node["node_attrs"]["clade_score"]["value"])
        if node["node_attrs"]["clade_score"]["value"]>60:
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
    with open(f'trial_{args.lineage}_clades.json', 'w') as fh:
        json.dump(data, fh)
