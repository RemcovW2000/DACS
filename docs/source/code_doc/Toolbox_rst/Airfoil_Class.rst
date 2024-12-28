Toolbox.Airfoil class
=====================

Main principle:
--------------
The airfoil class is built to do cross sectional analysis of an airfoil section
in a wing. Loads are assigned, and stresses are solved (shear and normal). During
this analysis the members in the airfoil section are discretized into boom and segment
objects, for which normal and shear stresses are calculated respectively.

During the failure analysis, for each member in the cross section, the worst
load case is considered, and using those loads, member failure analysis is carried
out.

Methods for calculating the shear flow can be found in the book 'Aircraft structures for engineering students',
in chapter 20. Authored by T.H.G. Megson. For this code, the structural idealisation was adapted to include
a large number of booms, which increases the accuracy of the shear flow calculation considerably.

Method docstrings:
------------------

.. automodule:: Toolbox.Airfoil
   :members:
   :undoc-members:
   :show-inheritance: