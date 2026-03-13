.. _sampling:

================
Sampling methods
================
This package works by calculating the rmsd against ground truth and applying an arbitrary cost function of self avoiding walks. A key part of this package rely on sampling self avoiding walks of a given length on a given lattice in a way that allows for observables to be estimated.

There are three sampling methods to choose from:
 * ``enumerate``: Iterate over all possible self avoiding walks of a given length ``n``. Note 	that the number of self avoiding walks is **exponential** in ``n``. Do not use it for large values of ``n`` (above ``10`` to ``20`` depending of the lattice and your computer).
 * ``dimerize``: The dimerize algorithm, performs an efficient uniformly random sampling of self avoiding walks.
 * ``pivot``: The pivot algorithm is a markovian process that samples self avoiding walks in an ergotic pseudo uniform way.


----------------
Pivot parameters
----------------
The pivot algorithm samples self avoiding walks (SAW) through an ergotic Markov chain. At the initialization of the algorithm and between each returned SAW, the algorithm generates and throw away several SAW. This process ensure that the SAW returned by the algorithm are independent from each others. Two parameters can be tuned:
 * ``thermalization_factor``: The pivot algorithm performs ``sequence_length * thermalization_factor`` accepted pivots before returning anything.
 * ``autocorrelation_factor``: The pivot algorithm performs ``autocorrelation * acceptation_rate`` steps (pivots accepted or rejected) between each returned SAW.

----------
References
----------
 1. Self avoinding walks: <https://en.wikipedia.org/wiki/Self-avoiding_walk>
 2. Methods to sample SAW: <https://arxiv.org/abs/hep-lat/9405016>
 3. Details on the Pivot algorithm: <https://clisby.net/projects/saw_feature/>
 4. Rust crate for the Pivot algorithm: <https://crates.io/crates/pivot_saw>