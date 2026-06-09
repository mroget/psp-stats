.. autofunction:: qpsp_stats.rmsd

>>> from qpsp_stats import rmsd
>>> sol = [[0.,0.,0.], [1.,1.,1.]]
>>> gt = [
... [[0.,0.,0.],[1.,-1.,1.]],
... [[1.,1.,1.],[0.,0.,0.]]
... ]
>>> rmsd(sol,gt)
4.458205583648681e-16