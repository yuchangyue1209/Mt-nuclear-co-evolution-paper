# Code provenance

The original copies of the 31 selected server scripts are preserved outside this public repository in the private project backup. Curated copies are placed under `scripts/` and use generic path placeholders. The private archive provides the audit trail between this repository and the working server while preventing personal paths and server details from entering the public release.

## Analysis mapping

| Repository directory | Analysis |
|---|---|
| `scripts/04_divergence_diversity` | codeml divergence and Pool-seq πN/πS analyses |
| `scripts/05_erc` | evolutionary-rate covariation |
| `scripts/06_pbs` | regional PBS and chromosome-XXI-controlled associations |
| `scripts/07_delta_af` | fixed shared-SNP ΔAF and recent-population evaluation |
| `scripts/08_geo_mtlineage` | geographic and mitochondrial-lineage models |
| `scripts/09_nonsynonymous` | amino-acid reconstruction, candidate prioritization, SIFT4G, and BLOSUM62 |
| `scripts/10_structure` | COX4I1 structural figure |

## Remaining audit

Curated scripts still contain paths from the original workstation and server. These paths must be replaced with configuration variables before public release. Numerical equivalence to the manuscript outputs must then be checked after each workflow is run.

## Neutral workflow labels

Public scripts use biological and analytical labels in place of personal names
in identifiers, path templates, comments, and filenames. Formal source citations
remain necessary: the direct/indirect/non-mitochondrial classification follows
Kuster et al. (2026), and gene-set ERC uses the gene-bootstrap framework of
Weaver et al. (2022). These citations do not label individual implementations.

Data and workspace roots remain configurable placeholders. Existing external
workspaces need not be renamed; configure the path variables to point to them.
Archived classification column headers are normalized on import by their suffix
and category values. Old ERC intermediate filenames and renamed plot inputs
must be regenerated using the consistently renamed upstream scripts, or staged
under the new filenames. Original working files are retained in the private backup.
