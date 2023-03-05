
import json, argparse

def add_provenance_clades(tree, key):
    '''
    Renames clades in the tree according to the clade_map.
    '''
    def rename_clade(n, parent_clade, current_clade, key):
        if key in n['branch_attrs'].get("labels",{}):
            if n["branch_attrs"]["labels"][key][0]!=current_clade:
                parent_clade = current_clade
                current_clade = n["branch_attrs"]["labels"][key][0]

            n["branch_attrs"]["labels"][key] = parent_clade + '-' + n["branch_attrs"]["labels"][key]

        n["node_attrs"][key] = {"value": parent_clade + '-' + n["node_attrs"][key]["value"]}
        if "children" in n:
            for c in n["children"]:
                rename_clade(c, parent_clade, current_clade, key)

    rename_clade(tree, 'none','none', key)


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Assign clades to a tree")
    parser.add_argument('--input', help='input josn')
    parser.add_argument('--key', default='new_clade', help="label to use for new clades")
    parser.add_argument('--output', default='renamed_tree.json')

    args = parser.parse_args()

    with open(args.input) as fh:
        data = json.load(fh)

    T = data["tree"]

    add_provenance_clades(T, args.key)

    with open(args.output, 'w') as fh:
        json.dump(data, fh)
