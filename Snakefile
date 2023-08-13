

rule fetch_and_rename:
    input:
        "config/{lineage}_old_to_new.json"
    output:
        "auspice/{prefix}renamed_{lineage}_{resolution}_{segment}.json"
    params:
        key = "renamed_clade"
    shell:
        """
        python3 scripts/rename_old_clades.py --lineage {wildcards.lineage} \
                --resolution {wildcards.resolution} \
                --segment {wildcards.segment} \
                --clade-map {input} \
                --key {params.key} \
                --output {output}
        """

rule suggest_new_clades:
    input:
        renamed_auspice = "auspice/{prefix}renamed_{lineage}_{resolution}_{segment}.json",
        clade_map = "config/{lineage}_old_to_new.json",
        weights = "config/weights.json",
    output:
        "auspice/{prefix}suggested_{lineage}_{resolution}_{segment}.json"
    params:
        old_key = "renamed_clade",
        new_key = "new_clade",
        add_to_existing = lambda w:"" if w.segment=='na' else "--add-to-existing"
    shell:
        """
        python3 scripts/add_new_clades.py --input {input.renamed_auspice} \
                --clade-map {input.clade_map} \
                --weights {input.weights} \
                --segment {wildcards.segment} \
                --lineage {wildcards.lineage} \
                {params.add_to_existing} \
                --old-key {params.old_key} --new-key {params.new_key} \
                --output {output}
        """


rule add_provenance:
    input:
        "auspice/{prefix}suggested_{lineage}_{resolution}_{segment}.json"
    output:
        "auspice/{prefix}with-provenance_{lineage}_{resolution}_{segment}.json"
    params:
        key = "new_clade"
    shell:
        """
        python3 scripts/add_provenance.py --input {input} --key {params.key}\
                --output {output}
        """

from datetime import date
today = date.today().strftime('%Y-%m-%d')
rule all:
    input:
        expand("auspice/{prefix}with-provenance_{lineage}_{resolution}_{segment}.json",
                prefix = [today+'_'],
                lineage = ['h3n2', 'h1n1pdm', 'vic'],
                resolution = ['2y', '6y'],
                segment = ['ha', 'na'])
