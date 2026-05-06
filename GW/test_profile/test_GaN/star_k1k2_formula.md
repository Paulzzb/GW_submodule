# Star(k1, k2) and Representative Notes

This file is intended to be edited frequently while refining the math definitions.

## 0. Symmetry Group and Action

Let

$$
\mathcal{G} = \{S_1, S_2, \dots, S_{N_s}\}
$$

be the set of symmetry matrices (point-group operations).

Define the action of $S \in \mathcal{G}$ on an ordered pair $(k_1, k_2)$ by

$$
S \cdot (k_1, k_2) := (S k_1, S k_2).
$$

## 1. Definition of Star(k1, k2)

For a given pair $(k_1, k_2)$, define

$$
\mathrm{Star}(k_1, k_2)
= \{(k_1', k_2') \mid \exists S \in \mathcal{G},\; k_1' = S k_1,\; k_2' = S k_2\}.
$$

This is exactly the orbit of $(k_1, k_2)$ under the group action of $\mathcal{G}$.

Recommended notation (abstract-algebra style):

$$
[(k_1, k_2)] := \mathrm{Star}(k_1, k_2).
$$

(If you prefer a shorter style in code comments, you can also write $[k_1\,k_2]$.)

## 2. Equivalence Relation

Define

$$
(k_1, k_2) \sim (\tilde{k}_1, \tilde{k}_2)
\iff \exists S \in \mathcal{G},\; (\tilde{k}_1, \tilde{k}_2) = (S k_1, S k_2).
$$

Pair-order convention used in this project:

$$
(k_1, k_2) \equiv (k_2, k_1),
$$

because the corresponding terms are elementwise conjugated. So pair order does not define a different class.

Then each class of $\sim$ is one orbit (one Star set):

$$
[(k_1, k_2)] = \{(\tilde{k}_1, \tilde{k}_2) : (\tilde{k}_1, \tilde{k}_2) \sim (k_1, k_2)\}.
$$

## 3. Representative of Each Orbit

For each orbit/class $[(k_1, k_2)]$, choose one representative

$$
\operatorname{rep}([(k_1, k_2)]) = (k_1^{\star}, k_2^{\star}) \in [(k_1, k_2)].
$$

Equivalently, define a representative set

$$
\mathcal{R} \subseteq \{(k_1, k_2)\}
$$

such that every orbit intersects $\mathcal{R}$ at exactly one point.

## 4. Representative Selection Rule (Editable)

Use this section to define your current practical rule.

Current rule (pseudo-code):

```matlab
marked = zeros(nk, nk);
rep_set = [];

for ik1 = 1:nk
  for ik2 = 1:nk
    if ~marked(ik1, ik2)
      rep_set = [rep_set; ik1, ik2];
      % Find all index pairs in Orbit(ik1, ik2) and mark them.
      % The orbit should include swapped order pairs as equivalent.
      % Prefer symmetry choices with G0 = 0 when possible.
    end
  end
end
```

## 5. Mapping Objects for Implementation

If $|\mathcal{R}| = N_{\mathrm{rep}}$, define:

- `k1k2_representation` of size `N_rep x 2`, storing representative pair indices.
- `k1k2_mapping(ik1, ik2, 1) = irep`, representative index.
- `k1k2_mapping(ik1, ik2, 2) = irot`, symmetry index satisfying
  $$
  k_1 = S_{irot} k_{1,\mathrm{rep}} + G_0, \quad
  k_2 = S_{irot} k_{2,\mathrm{rep}} + G_0.
  $$

Here $G_0$ is forced to be the same reciprocal-lattice shift for both components.

Practical preference for choosing $S_{irot}$:

1. If possible, choose $S_{irot}$ such that $G_0 = 0$.
2. Otherwise choose any valid $S_{irot}$ with a shared nonzero $G_0$.

## 6. TODO

- [x] Pair order does not matter: $(k_1,k_2) \equiv (k_2,k_1)$.
- [x] Same reciprocal shift is enforced: $G_1=G_2=G_0$.
- [x] Finalize orbit-marking implementation details in MATLAB.
- [x] Sync notation with MATLAB variable names in demo scripts.
