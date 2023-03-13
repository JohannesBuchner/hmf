=======================================
Halo Mass Function in Python
=======================================

One would think that in 2018 there would be a python package that computes 
the halo mass function at an arbitrary redshift.
Well, there isn't.

This is a stupid replacement that plots them from large N-body simulations using a lookup table.

See code for documentation and cite:

https://www.cosmosim.org/cms/simulations/smdpl/
https://www.cosmosim.org/cms/simulations/mdpl2/
https://www.cosmosim.org/cms/simulations/bigmdpl/

Example outputs
------------------------

Halo mass function:

.. image:: https://raw.githubusercontent.com/JohannesBuchner/hmf/master/hmf_M.png

Cumulative halo mass function:

.. image:: https://raw.githubusercontent.com/JohannesBuchner/hmf/master/hmf_cumulative_M.png

Redshift evolution above some mass cut:

.. image:: https://raw.githubusercontent.com/JohannesBuchner/hmf/master/hmf_cumulative_z.png

See also
-----------

* http://roban.github.com/CosmoloPy/
* http://halotools.readthedocs.io/
* http://yt-project.org/  (removed its halo mass function module because of code rot)
* http://docs.astropy.org/en/stable/cosmology/
* https://hmf.readthedocs.io/  (computes older approximations, many parameters)
