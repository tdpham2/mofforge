"""Check periodic clustering against exhaustive selection and report timings."""

from time import perf_counter

import numpy as np
from pymatgen.core import Lattice

from mofforge.adsorbate.sites import _cluster_void_points


def main():
    rng = np.random.default_rng(31)
    lattice = Lattice.from_parameters(20, 22, 24, 75, 80, 65)
    for n in (250, 1000, 4000):
        frac = rng.random((n, 3))
        scores = rng.random(n)
        started = perf_counter()
        chosen = []
        for i in np.argsort(-scores, kind="stable"):
            if not chosen or lattice.get_all_distances([frac[i]], frac[chosen]).min() >= 1:
                chosen.append(i)
        exhaustive_time = perf_counter() - started
        started = perf_counter()
        actual = _cluster_void_points(frac, lattice.get_cartesian_coords(frac), scores, lattice, 1)
        indexed_time = perf_counter() - started
        np.testing.assert_allclose([s.frac_coords for s in actual], frac[chosen])
        print(
            f"{n:5d} points: exhaustive={exhaustive_time:.3f}s indexed={indexed_time:.3f}s "
            f"centers={len(actual)}"
        )


if __name__ == "__main__":
    main()
