Docker
======

.. rst-class:: lead

   Run NNV inside a Docker container with MATLAB.

----

Requirements
------------

- ~20 GB disk space (18.2 GB image)
- Valid MATLAB license

Build the Image
---------------

Two variables can be configured in the ``Dockerfile``:

1. **MATLAB release** (default R2024b):

   .. code-block:: dockerfile

      ARG MATLAB_RELEASE=R2024b

2. **License server** (port@hostname format):

   .. code-block:: dockerfile

      ARG LICENSE_SERVER="27009@licenseserver.it.vanderbilt.edu"

Then build:

.. code-block:: bash

   docker build . -t nnv

Run Interactively
-----------------

.. code-block:: bash

   docker run -it nnv

This opens an interactive MATLAB session with NNV installed and ready to use.
