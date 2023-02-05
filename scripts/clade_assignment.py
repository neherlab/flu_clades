import json, argparse
import urllib.request, json
from collections import defaultdict
import numpy as np
from config.weights import weights
from config.clade_map import old_to_new_clades

def get_existing_clade_labels(tree):
    '''
    Returns a list of existing clade labels in the tree.
    '''
    def add_clades(n, clades):
        if "labels" in n["branch_attrs"] and "clade" in n["branch_attrs"]["labels"]:
            clades.add(n["branch_attrs"]["labels"]["clade"])
        if "children" in n:
            for c in n["children"]:
                add_clades(c, clades)

    clades = set()
    add_clades(tree, clades)
    return clades


def rename_clades(tree, clade_map):
    '''
    Renames clades in the tree according to the clade_map.
    '''
    def rename_clade(n, clade_map):
        if "labels" in n["branch_attrs"] and "clade" in n["branch_attrs"]["labels"]:
            clade = n["branch_attrs"]["labels"]["clade"]
            if clade in clade_map:
                n["branch_attrs"]["labels"]["clade"] = clade_map[clade][0]
                n["node_attrs"]["clade"] = {"value":clade_map[clade][1]}

        clade = n["node_attrs"]["clade"]["value"]

        if "children" in n:
            for c in n["children"]:
                c["node_attrs"]["clade"] = {"value": clade}
                rename_clade(c, clade_map)

    rename_clade(tree, clade_map)

def label_backbone(tree):
    def label_backbone_recursive(n):
        on_backbone = False
        if "children" in n:
            for c in n["children"]:
                label_backbone_recursive(c)
                if c["backbone"]:
                    on_backbone = True

        if on_backbone or ("labels" in n["branch_attrs"] and "clade" in n["branch_attrs"]["labels"]):
            n["backbone"] = True
        else:
            n["backbone"] = False

    label_backbone_recursive(tree)

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

def gather_tip_counts(n, attenuate=None, max_value=0, ignore_backbone=False):
    '''
    Assigns a node attribute "ntips" to each node that counts the number of downstream tips.
    This can be attenuated similar to the LBI calculation.
    '''
    if "children" in n:
        tmp = 0
        for c in n["children"]:
            tmp_max = gather_tip_counts(c, attenuate=attenuate, max_value=max_value, ignore_backbone=ignore_backbone)
            if (ignore_backbone and c["backbone"])==False:
                tmp += c["ntips"]*attenuate(c)
            else:
                print("Ignoring clade breakpoint at node", c['backbone'])
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

def rename_new_clades(n, name_map):
    n['node_attrs']['trial_clade']['value'] = name_map(n['node_attrs']['trial_clade']['value'])
    if 'children' in n:
        for c in n['children']:
            rename_new_clades(c, name_map)

def rename_clade_labels(n, name_map):
    if "labels" in n['branch_attrs'] and 'trial_clade' in n['branch_attrs']['labels']:
        print(n['branch_attrs']['labels']['trial_clade'])
        n['branch_attrs']['labels']['trial_clade'] = name_map(n['branch_attrs']['labels']['trial_clade'])
    if 'children' in n:
        for c in n['children']:
            rename_clade_labels(c, name_map)

def score(n, weights=None, ntip_scale=1, ignore_backbone=False):
    if weights is None:
        weights = {}

    if ignore_backbone and n['backbone']:
        return 0.0
    score = np.sqrt(n['ntips']/ntip_scale)
    mut_weight = sum([weights.get(int(x[1:-1]),1) for x in n['branch_attrs']['mutations'].get('HA1',[])])
    if mut_weight>0 and ('clade' not in n['branch_attrs'].get('labels', {})):
        score += mut_weight/(4 + mut_weight)
    else:
        score = 0.0
    return score


def assign_score(n, score=None, weights=None, ntip_scale=1, ignore_backbone=False):
    if "children" in n:
        for c in n["children"]:
            assign_score(c, score=score, weights=weights, ntip_scale=ntip_scale, ignore_backbone=ignore_backbone)
    n['node_attrs']["score"] = {'value': score(n, weights=weights, ntip_scale=ntip_scale, ignore_backbone=ignore_backbone)}


def assign_clade(n, clade):
    if "children" in n:
        for c in n["children"]:
            assign_clade(c, clade)
    n['node_attrs']["clade"] = {'value': clade}

def assign_clades_to_branch(n, hierarchy, cutoff=1.0):
    node = first_above(n, 'score', cutoff)
    if node is not None:
        if "clade_score" not in node["node_attrs"]:
            print("new clade",node['name'], node['div'], node['ntips'], node['node_attrs']['score']['value'])

            if 'labels' not in node['branch_attrs']:
                node['branch_attrs']['labels'] = {}

            node["node_attrs"]["clade_score"] = {"value":node["node_attrs"]["score"]["value"]}
            print(node["node_attrs"]["clade_score"])

            parent_clade = node['node_attrs']["clade"]["value"]
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
    T["node_attrs"]["clade"] = {"value": "none"}
    rename_clades(T, old_to_new_clades)
    existing_clades = get_existing_clade_labels(T)
    label_backbone(T)
    max_value = gather_tip_counts(T, lambda x:np.exp(-len(x['branch_attrs']['mutations'].get('HA1',[]))/1),
                                  ignore_backbone=True)

    root_clade = ('A',)
    T['div']=0
    assign_divergence(T)
    assign_score(T, score, weights=weights[args.lineage], ntip_scale=max_value, ignore_backbone=True)

    assign_clade(T, root_clade)
    hierarchy = defaultdict(list)
    for clade_name, full_clade in old_to_new_clades.values():
        if len(full_clade)==1:
            continue
        hierarchy[full_clade[:-1]].append(full_clade[-1])
        hierarchy[full_clade] = []

    hierarchy = assign_clades_to_branch(T, hierarchy, cutoff=0.75)

    aliases = {root_clade: "A"}
    major_clades = sorted(set([x[0] for x in existing_clades]))
    alpha = 'ABCDEFGHIJKLMNOPQRST'
    for clade in clades:
        node  = clades[clade]
        print(node["node_attrs"]["clade_score"]["value"])
        if node["node_attrs"]["clade_score"]["value"]>1:
            year = int(node["node_attrs"]["num_date"]["value"]) - 2000
            if len(major_clades)>len(alpha):
                letter = str(len(major_clades)//len(alpha))
            else: letter=''
            letter += alpha[len(major_clades)%len(alpha)]
            aliases[node["node_attrs"]["trial_clade"]["value"]] = f"{letter}"
            major_clades.append(letter)


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

    rename_new_clades(T, lambda x:name_map[x])
    rename_clade_labels(T, lambda x:name_map[x])

    data['meta']['colorings'].append({'key':"trial_clade", 'type':'ordinal', 'title':"trial clade"})
    data['meta']['colorings'].append({'key':"score", 'type':'continuous', 'title':"clade score"})
    with open(f'auspice/trial_{args.lineage}_clades.json', 'w') as fh:
        json.dump(data, fh)
