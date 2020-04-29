Welcome to Image3d's documentation!
========================================

Image3d is a python tools created to work on 3D images organized on regular voxel such as the one obtained in tomography. It is an open source code under CC-BY-CC licence (https://creativecommons.org/licenses/by/2.0/fr/deed.en). There is no guarantee if you are using it. It has been tested with python 3.7


.. toctree::
   :maxdepth: 1
   :caption: Documentation:


Installation
============

For simple user
***************

If you want just to use the toolbox without doing any devellopement

You need **setvector3d** to use **image3d**

.. code:: bash

	pip install git+https://github.com/ThomasChauve/setvector3d
	pip install git+https://github.com/ThomasChauve/image3d


For develloper
**************

If you want to access the code and be able to do some modification

.. code:: bash

    git clone https://github.com/ThomasChauve/setvector3d
    cd setvector3d/
    pip install .

.. code:: bash

    git clone https://github.com/ThomasChauve/image3d
    cd image3d/
    pip install .


Then you will find all the package in python using

.. code:: python

    import setvector3d
    import image3d

Uninstall
*********

.. code:: bash
    
    pip uninstall setvector3d
    pip uninstall image3d


.. toctree::
    :maxdepth: 1
    :numbered:
    :caption: Documentation

    Doc1/Exemples
    Doc2/Devellopment

.. toctree::
    :maxdepth: 1
    :numbered:
    :caption: CLASS

    object


Contact
=======
:Author: Thomas Chauve
:Contact: thomas.chauve@univ-grenoble-alpes.fr

:organization: UiO
:status: This is a "work in progress"
:version:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
