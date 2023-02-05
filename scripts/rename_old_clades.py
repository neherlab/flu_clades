
import json, argparse
import urllib.request, json

def rename_clades(tree, clade_map, key):
    '''
    Renames clades in the tree according to the clade_map.
    '''
    def rename_clade(n, clade_map, key):
        if "labels" in n["branch_attrs"] and "clade" in n["branch_attrs"]["labels"]:
            clade = n["branch_attrs"]["labels"]["clade"]
            if clade in clade_map:
                n["branch_attrs"]["labels"][key] = clade_map[clade][0]
                n["node_attrs"][key] = {"value":clade_map[clade][0]}
                n["node_attrs"][f"full_{key}"] = {"value":clade_map[clade][1]}

        clade = n["node_attrs"][key]["value"]

        if "children" in n:
            for c in n["children"]:
                c["node_attrs"][key] = n["node_attrs"][key]
                c["node_attrs"][f"full_{key}"] = n["node_attrs"][f"full_{key}"]

                rename_clade(c, clade_map, key)

    tree["node_attrs"][key] = {"value":"none"}
    tree["node_attrs"][f"full_{key}"] = {"value":tuple()}
    rename_clade(tree, clade_map, key)


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Assign clades to a tree")
    parser.add_argument('--lineage', default='h3n2', help="lineage to assign clades to")
    parser.add_argument('--clade-map',help="json with renamed clades")
    parser.add_argument('--key', default='new_clade', help="label to use for new clades")
    parser.add_argument('--resolution', default='6y')
    parser.add_argument('--output', default='renamed_tree.json')

    args = parser.parse_args()
    with open(args.clade_map) as fh:
        old_to_new_clades = json.load(fh)

    with urllib.request.urlopen(f"https://nextstrain.org/charon/getDataset?prefix=flu/seasonal/{args.lineage}/ha/{args.resolution}") as url:
        data = json.loads(url.read())

    T = data["tree"]

    rename_clades(T, old_to_new_clades, args.key)

    data['meta']['colorings'].append({'key':args.key, 'type':'ordinal', 'title':args.key})

    with open(args.output, 'w') as fh:
        json.dump(data, fh)
