
import json, argparse
import urllib.request, json

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Assign clades to a tree")
    parser.add_argument('--lineage', default='h3n2', help="lineage to assign clades to")
    parser.add_argument('--clade-map',help="json with renamed clades")
    parser.add_argument('--output', default='renamed_clades_tsv.json')

    args = parser.parse_args()
    with open(args.clade_map) as fh:
        old_to_new_clades = json.load(fh)

    suffix = ""
    with urllib.request.urlopen(f"https://raw.githubusercontent.com/nextstrain/seasonal-flu/master/config/{args.lineage}/ha/clades{suffix}.tsv") as url:
        data = url.read().decode().split('\n')

    with open(args.output, "w") as fh:
        for line in data:
            entries = line.split("\t")
            if len(entries) >= 3:
                for i in [0,2]:
                    if entries[i] in old_to_new_clades:
                            entries[i] = old_to_new_clades[entries[i]][0]

            fh.write('\t'.join(entries)+'\n')

