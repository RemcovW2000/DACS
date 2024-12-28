Toolbox.Member class
====================

Main principle:
--------------

A member is a panel with a certain size (a and b) where a is the width, and b is the depth of the member.
The panel is either a laminate or a sandwich at this moment.

Whereas a panel object can only fail due to first ply failure (laminate and sandwich panels) or wrinkling and
dimpling (only sandwich panels), a member has an emergent failure mode, namely buckling.

Loads are assigned to the member, optionally the curvature of the member is assigned, the critical load
intensities are calculated at which buckling occurs.

Buckling equations are obtained from the book 'Design and analysis of composite structures' authored by
Dr. Christos Kassapoglou.

Method docstrings:
------------------

.. automodule:: Toolbox.Member
   :members:
   :undoc-members:
   :show-inheritance:
