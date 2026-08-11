# Validation Plan for Symmetry/Product Consistency

## Scope

This document defines two validation targets for wave-function symmetry and pair-product equivariance.

## Constraint

- Do not create a new SAVE directory.
- Reuse existing data and existing stage/snapshot artifacts under test_Si.

## Validation 1: Fixed-k Symmetry Mapping

### Statement (Fixed-k)

If S k' = k (including k = k'), then wave functions must satisfy the corresponding symmetry relation.

### Key example

For k = 0 and k' = 0, any symmetry rotation S that belongs to the crystal symmetry group should preserve the wave-function structure at Gamma, i.e. the wave function has rotational symmetry under all valid S.

### Check idea (Fixed-k)

1. Enumerate symmetry rotations S from restored symmetry manager data.
2. For each S and each selected band/spin channel at k = 0:
   - rotate psi(k=0) in real/reciprocal representation,
   - compare with the expected psi at mapped k' (still k'=0 for Gamma case),
   - evaluate residual norm and phase-aligned mismatch.
3. Report max/mean residual and pass-fail against tolerance.

## Validation 2: Product-Rotation Commutativity on Equivalent Pairs

### Statement (Pair-product)

For equivalent pairs (k1', k2') ~ (k1, k2), with k1' = S k1 and k2' = S k2:

- multiply first then rotate:
  rotate( psi(k1) * psi(k2) )

must equal

- rotate first then multiply:
  psi(k1') * psi(k2')

within numerical tolerance.

### Check idea (Pair-product)

1. Use existing k-pair mapping data to identify representative and mapped pairs.
2. For each tested pair and rotation S:
   - construct product field from original pair then apply S,
   - construct product field from rotated/mapped pair directly,
   - compare residual norm after optional phase alignment.
3. Aggregate statistics (max/mean residual, fail count).

## Expected Deliverables (next steps)

1. Demo entry script for these validations.
2. Supporting helper functions.
3. Runtime log file summarizing pass/fail and residual statistics.
