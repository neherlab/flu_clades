
import json, argparse
import urllib.request, json
from rename_old_clades import rename_clades


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Assign clades to a tree")
    parser.add_argument('--lineage', default='rsv/a', help="lineage to assign clades to")
    parser.add_argument('--clade-map',help="json with renamed clades")
    parser.add_argument('--key', default='new_clade', help="label to use for new clades")
    parser.add_argument('--output', default='renamed_tree.json')

    args = parser.parse_args()
    with open(args.clade_map) as fh:
        old_to_new_clades = json.load(fh)

    with urllib.request.urlopen(f"https://nextstrain.org/charon/getDataset?prefix={args.lineage.replace('-','/')}/genome/") as url:
        data = json.loads(url.read())

    T = data["tree"]

    rename_clades(T, old_to_new_clades, args.key)

    data['meta']['colorings'].append({'key':args.key, 'type':'ordinal', 'title':args.key})

    with open(args.output, 'w') as fh:
        json.dump(data, fh)
