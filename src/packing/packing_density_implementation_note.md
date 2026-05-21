# Packing Density Implementation Note

## Current implementation

The current `src/packing/density.py` implementation calculates packing density in native Python.

For each selected BDamage atom, the code loops over every atom in the trimmed crystal neighbour cloud and counts how many neighbours fall within the packing-density threshold.

Conceptually:

```text
for each selected BDamage atom:
    count = 0
    for each trimmed neighbour atom:
        calculate distance between selected atom and neighbour atom
        if distance < packing_density_threshold:
            count += 1
```

The implementation compares squared distances rather than true distances:

```text
dx = selected_x - neighbour_x
dy = selected_y - neighbour_y
dz = selected_z - neighbour_z

distance_squared = dx*dx + dy*dy + dz*dz
```

Then it checks:

```text
distance_squared < packing_density_threshold * packing_density_threshold
```

This avoids unnecessary square-root calculations while giving the same inclusion result as comparing actual Euclidean distances.
The final packing-density value subtracts one from this raw count to match
RABDAM 2 behaviour, where the selected atom's own central-cell copy is assumed
to be present in the trimmed neighbour cloud and is removed from the contact
count.

## Why this is acceptable for the baseline

The native-Python implementation is simple, explicit, and easy to test. It makes the packing-density definition transparent:

```text
packing density = number of trimmed crystal atoms within the threshold distance
                  of a selected asymmetric-unit atom
```

This is useful during the first RABDAM 3 BDamage implementation because correctness and compatibility are more important than speed.

The current code deliberately avoids more complex neighbour-search optimizations so that the calculation is easy to reason about while validating RABDAM 3 against previous RABDAM behavior.

## Performance limitation

The current algorithm is an all-against-all search between selected atoms and trimmed neighbour atoms.

If there are:

```text
N selected atoms
M trimmed neighbour atoms
```

then the calculation performs approximately:

```text
N * M distance checks
```

This may become slow for large structures or large neighbour clouds.

## Possible future optimizations

After the baseline BDamage implementation is validated, packing-density counting could be optimized using one of the following approaches:

- NumPy vectorized distance calculations
- `scipy.spatial.cKDTree`
- cell lists / spatial hashing
- another spatial-neighbour search structure

Any optimization should preserve the same inclusion rule:

```text
distance < packing_density_threshold
```

and should be tested against the native-Python implementation to confirm that neighbour counts remain identical.

## Recommended future policy

Keep the native-Python implementation as the reference implementation, even if a faster backend is added later.

A future optimized implementation could be added behind an option such as:

```text
method = "python" | "kdtree" | "numpy"
```

The native-Python method would remain useful for testing, debugging, and confirming numerical compatibility.

## Summary

The current packing-density calculation is intentionally simple:

```text
selected atoms + trimmed neighbour cloud + threshold
→ count neighbours within threshold for each selected atom
```

It is not the fastest possible approach, but it is transparent, testable, and appropriate for the first compatibility-focused RABDAM 3 BDamage implementation.
