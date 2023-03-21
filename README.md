# Pipeline for automated suggestion of new clades of seasonal influenza viruses


The general idea of this workflow is to fetch a current tree from Nextstrain that is already annotated with currently used clades, search for clusters in the tree that meet certain criteria (see below), and suggest the branches leading to these clusters as novel clades.

## Clade naming criteria

The purpose of this workflow is to be run at regular intervals and suggest new clades to a tree that already has annotation.
It is thus not meant to find a globally optimal partition of an unlabeled tree, but to identify novel groups that meet designation criteria as they grow and evolve over time.
The three basic aspects that enter the criterion are
 * **size:** big groups should have a higher priority for designation.
 * **divergence:** the more mutations have accumulated relative to the break point of the parent clade, the higher the priority of a novel clade
 * **specific mutations:** Ideally, breakpoints sit on long branches with significant mutation.

There are many different ways in which these aspects could be combined into a single score.
The proposal below seems to work well.

### Size and phylogenetic structure

To single out clades in the tree that are of sufficient size and are phylogenetically distinct, we propose a measure similar to the LBI.
$$
\phi_n = \sum_{c\in n} d(1-e^{l_c/d}) +  \phi_c e^{l_c/d}
$$
where $l_c$ is the length of the branch leading to child node $c$ and $d$ is a distance scale determining how rapidly this score is "forgotten" along the tree.
Terminal nodes have $\phi_c=0$.
To make the values of the score independent of the size of the tree, we normalize the score to its maximal value and take the square root.
$$
\psi_n = \sqrt{\phi_n/\phi_{max}}
$$
This way, the value of this score varies between 0 and 1 with a clade that a quarter of the maximal clade being assigned a value of 0.5.

### Branch value.
To prioritize branches with important mutations, we assign each mutation a weight and sum these weights for all mutations on the branch.
For A/H3N2 HA, this is currently 3 for each of the Koel-sites, 2 for epitope sites, 1 for other HA1 or HA2 mutations, and 0 for synonymous mutations.
These weights can be specified in a position specific way for each lineage and segment.
The sum of weights $w_n$ of a branch $n$ is transformed to a value between 0 and 1 via
$$
\beta_n = \frac{w_n}{w_n+4}
$$


### Divergence

Divergence is a simple count of amino acid mutations in a subset of proteins (HA1 for influenza) along the tree since the last clade breakpoint.
For each node $n$, this count $b_n$ is transformed to a value between 0 and 1 via
$$
\delta_n = \psi_n\frac{b_n}{4+b_n}
$$
The multiplication of the divergence with the phylogenetic score ensure that assignments triggered by the divergence criteria tend to mark sizeable clades.
Otherwise, gradual increase of divergence would trigger clade designations in undesired places.

## Clade assignment

New clades are assigned by walking through the tree in pre-order (parents before children) and a new clade is designated if the sum of scores exceeds a certain threshold $\theta$.
$$
\psi_n + \beta_n + \delta_b > \theta_n
$$
The current value for the threshold in the A/H3N2 clade assignment is $\theta=0.85$.
This means that mutations on a branch alone are insufficient to trigger a clade assignment unless their weight sums to more than 22.
But together with a phylogenetic score, the threshold can be crossed.

