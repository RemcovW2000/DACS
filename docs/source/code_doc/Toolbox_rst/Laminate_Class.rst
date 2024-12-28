Toolbox.Laminate class
=====================

Main principle:
--------------
A laminate is build up out of individual lamina. These lamina can individually fail. When one lamina
fails, this is referred to as 'first ply failure'. This class serves to calculate the load intensities
of individual lamina in a laminate, given the load intensities applied on the total laminate. Then,
a failure analysis can be done for each individual lamina. The critical load intensity at which first ply
failure occurs can be calculated.

A Laminate object can be used by the Member class as the 'panel' attribute.

Method docstrings:
------------------
.. automodule:: Toolbox.Laminate
   :members:
   :undoc-members:
   :show-inheritance:
