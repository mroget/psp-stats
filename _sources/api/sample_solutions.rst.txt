===========================
qpsp_stats.sample_solutions
===========================

.. autofunction:: qpsp_stats.sample_solutions


"""""""""""""""
Sampling method
"""""""""""""""
``method`` can take the values:
 * ``"pivot"``
 * ``"dimerize"``
 * ``"enumerate"``

See section :ref:`sampling` for a detailed description of these methods.


"""""""
Lattice
"""""""
``lattice`` can take the values:
 * ``"tetrahedral"``: The thetrahedral lattice (degree 4).
 * ``"cubic"``: The cubic lattice (degree 6).
 * ``"bcc"``: The Body Centered Cubic (BCC) lattice (degree 8).
 * ``"fcc"``: The Face Centered Cubic (FCC) lattice (degree 12).