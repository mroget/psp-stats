.. _cost:

==============
Cost functions
==============
This package works by calculating the rmsd against ground truth and applying an arbitrary cost function of self avoiding walks. In particular, there is three ways to pass an arbitrary cost function to this package: 
 * A Python function
 * A Dictionnary
 * A string tuple containing lua code

---------------------
Use a python function
---------------------
This is the most convenient way to pass an arbitrary cost function. It is also the slowest since it executes python code directly.

The python function you put in arguments must have the following signature:
``fun(sequence : str, structure : list) -> float``,
where ``structure`` is of shape ``(len(sequence), 3)`` elements of type ``float``.


-----------------
Use a dictionnary
-----------------
This method provides the fastest way to calculate statistics, but is also quite restrictive in which cost function can be passed in arguments.

This method assume that the function ``cost`` that you want to pass has the following formula:

.. math::

	cost(\text{seq}, \text{structure}, w) = \sum_{0 \leq i < j < \text{len}(\text{seq})} \epsilon_{i,j}(\text{seq},\text{structure},w),

where 

.. math::

	\epsilon_{i,j}(\text{seq},\text{structure},w) = \left\{
	\begin{matrix}
		w_{\text{seq}_i,\text{seq}_j} & \text{if } j-i \geq k_{\text{min}} \text{ and } \Vert \text{structure}_i - \text{structure}_j \Vert_2 \leq d_\text{max}\\
		0 & \text{otherwise}\\
	\end{matrix}
		\right. .


This is essentially a sum of all the contact energies. All possible contact energies are described in dictionnary ``w`` which must be passed as the cost argument. Contacts are accounted for in the total cost function only if the euclidian distance is close enough and if the sequence distance is large enough. ``kmin`` and ``dmax`` are optional arguments.

------------
Use lua code
------------
This method allows to pass any cost function coded in lua. It is also much more efficient than passing a python function since the lua code is compiled JIT. This allows perfomances close to compiled code. 

To use this method, you must give a tuple of strings ``(lua_code, function_name)``. ``lua_code`` is a string containing the lua code. You can add as many declarations as you want, but the code must contain a function named ``function_name``. The cost will be calculated by calling ``function_name`` on couples of sequence and structure. Thus the function ``function_name`` must have the signature 
``fun(sequence : str, structure : list) -> float``,
where ``structure`` is of shape ``(len(sequence), 3)`` elements of type ``float``.

------
Remark
------
If you pass some python or some lua code, this code will be executed as is. If this code returns an error, the package will forward the error and stop.