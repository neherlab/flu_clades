# Pipeline for automated suggestion of new clades of seasonal influenza viruses


The general idea of this workflow is to fetch a current tree from Nextstrain that is already annotated with currently used clades, search for clusters in the tree that meet certain criteria (see below), and suggest the branches leading to these clusters as novel clades.

## Clade naming criteria

The purpose of this workflow is to be run at regular intervals and suggest new clades to a tree that already has annotation.
It is thus not meant to find a globally optimal partition of an unlabeled tree, but to identify novel groups that meet designation criteria as they grow and evolve over time.
The three basic aspects that enter the criterion are
 * **size:** big groups should have a higher priority for designation.
 * **divergence:** the more mutations have accumulated relative to the break point of the parent clade, the higher the priority of a novel clade
 * **specific mutations:** Ideally, breakpoints sit on long branches with significant mutation. Such mutations will be better defined for well studied segments/genomes like the HA segment of H3N2 influenza viruses.

There are many different ways in which these aspects could be combined into a single score.
The proposal below seems to work well.

### Size and phylogenetic structure

To single out clades in the tree that are of sufficient size and are phylogenetically distinct, we propose a measure similar to the LBI.
```math
\phi_n = \sum_{c\in n} (1-e^{-l_c/d}) +  \phi_c e^{-l_c/d}
```
where $l_c$ is the length of the branch leading to child node $c$ and $d$ is a distance scale determining how rapidly this score is "forgotten" along the tree.
Terminal nodes have $\phi_c=1$.
This way, the contribution of each terminal to its parent $(1-e^{-l_c/d}) + e^{-l_c/d} = 1$ regardless of its branch length.
Each internal branch has a maximum contribution of 1, less if the branch is short.
In addition, we calculate the number of leafs $l_n$ below each branch in the tree.

### Branch score
To prioritize branches with important mutations and to pick branches with many amino acid substitutions over those with only synonymous changes, we assign each mutation a weight and sum these weights for all mutations on the branch.
For A/H3N2 HA, this is currently 3 for each of the Koel-sites, 2 for epitope sites, 1 for other HA1 or HA2 mutations, and 0 for synonymous mutations.
These weights can be specified in a position specific way for each lineage and segment.
The branch score is called $b_n$ below.

### Divergence

Divergence is a simple count of amino acid mutations in a subset of proteins (HA1 for influenza) along the tree since the last clade breakpoint.
The more divergence has accumulated, since the last break point, the more readily a new clade should be designated.
The divergence score is called $d_n$ below.

## Aggregation of scores
Above, we defined three scores that together should determine whether a branch should be designated the root of a novel clade of lineage.
The branch score and the divergence score essentially count mutations and increase with the depth of the tree, but are fairly independent of sampling density.
The phylogenetic bushiness $\phi_c$ and the tip count depend strongly on sampling density.
To consistently combine these into a binary outcome, we need to normalize them.
The suggested normalization for the divergence, the branch score, and the bushiness score is linear saturation.
If $x$ is the raw score, we calculate the normalized score $X$ as
```math
X = \frac{x}{x+s_x}
```
where $s_x$ is a scale associated with the score.
For the branch score and the divergence, this is essentially given by the number of mutations considered a significant.
Since the phylogenetic score depends on the sample size, it makes sense to choose that normalization to be proportional to the tree, which is conveniently done by picking the max or median value.
We currently pick the median.

The three normalized scores are then added and if their sum exceeds a cut-off, a new clade is added.
```math
\Phi_n + B_n + D_n > C
```
Not that the divergence of the downstream part of the tree need recalculating after designating a new clade.
The clade size criterion is added as a hard binary threshold.

## Clade assignment

New clades are assigned by walking through the tree in pre-order (parents before children) and a new clade is designated if the sum of scores exceeds a certain threshold $\theta$.
```math
\psi_n + \beta_n + \delta_b > \theta
```
The current value for the threshold in the A/H3N2 clade assignment is $\theta=1.0$.
This means any one score alone is insufficient to trigger a new lineage.
But together with a phylogenetic score and the divergence contribution, the threshold can be crossed for a single HA1 mutation.

