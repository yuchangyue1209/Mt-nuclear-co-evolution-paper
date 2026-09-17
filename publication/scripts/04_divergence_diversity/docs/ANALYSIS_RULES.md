# Analysis rules

## Zero rates

| Observed state | dN | dS | omega |
|---|---:|---:|---:|
| invariant | 0 | 0 | NA |
| synonymous only | 0 | >0 | 0 |
| nonsynonymous only | >0 | 0 | NA |
| both | >0 | >0 | estimated |

Never add a pseudocount to dS. PAML boundary values `0.0001` and `999` are flags, not literal biological estimates.

## Primary analyses

- dN and dS: all 17,965 quality-eligible genes, including invariant zeros.
- omega ES1: genes with at least one expected synonymous substitution and no upper-bound estimate.
- omega ES2: stricter sensitivity analysis requiring at least two expected synonymous substitutions.
- Use dN as the primary response, dS as a mutation-rate/denominator check, and omega as supporting evidence.

## Current biological result

- Direct n-mt genes are not uniformly accelerated.
- Nmt-RP is faster than cyto-RP but not faster than genome-wide background.
- Nmt-ARS is faster than cyto-ARS and genome-wide background.
- Nuclear core subunits show a preliminary elevation relative to noncore subunits.

