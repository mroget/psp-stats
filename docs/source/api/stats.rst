================
qpsp_stats.stats
================

.. autofunction:: qpsp_stats.stats



"""""""""""""""
Sampling method
"""""""""""""""
``method`` can take the values:
 * ``"pivot"``
 * ``"dimerize"``
 * ``"iterate"``

See section :ref:`sampling` for a detailed description of these methods.


"""""""
Lattice
"""""""
``lattice`` can take the values:
 * ``"tetrahedral"``: The thetrahedral lattice (degree 4).
 * ``"cubic"``: The cubic lattice (degree 6).
 * ``"bcc"``: The Body Centered Cubic (BCC) lattice (degree 8).
 * ``"fcc"``: The Face Centered Cubic (FCC) lattice (degree 12).


"""""""""""""
Cost function
"""""""""""""
``cost`` can be of three different types:
 * A python function: Definitely the most convenient way to pass a cost function but slow since python code is executed.
 * A dictionnary ``d``: For any pair of (a1, a2) of amino acids, ``d[a1][a2]`` must be a float reprsenting the energy contribution of this pair. Fast but restrictive.
 * A tuple ``(lua_code, function_name)``: ``lua_code`` must contain a function named ``function_name``. Executed with luajit which is much faster than python.

See section :ref:`cost` for more information on how to pass a cost function.