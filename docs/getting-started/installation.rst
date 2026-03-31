Installation
============

.. rst-class:: lead

   NNV requires MATLAB 2023a or newer with several toolboxes. This guide
   covers installation on all platforms.

----

Prerequisites
-------------

**MATLAB 2023a or newer** with the following toolboxes:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Toolbox
     - Used For
   * - Computer Vision
     - Image processing operations
   * - Control Systems
     - Control system plant models
   * - Deep Learning
     - Neural network layers and operations
   * - Image Processing
     - Image manipulation and analysis
   * - Optimization
     - LP solving (linprog)
   * - Parallel Computing
     - Multi-core reachability analysis
   * - Statistics and Machine Learning
     - Statistical operations
   * - Symbolic Math
     - Symbolic computation
   * - System Identification
     - System modeling utilities

Clone the Repository
--------------------

NNV depends on several submodules (CORA, NNMT, HyST, onnx2nnv), so you must
clone recursively:

.. code-block:: bash

   git clone --recursive https://github.com/verivital/nnv.git

Ubuntu Quick Install
--------------------

On Ubuntu, a one-step installation script is provided:

.. code-block:: bash

   chmod +x install_ubuntu.sh
   ./install_ubuntu.sh

Then skip to :ref:`verify-installation`.

Windows / macOS Installation
----------------------------

1. Install MATLAB (2023a or newer) with the toolboxes listed above.

2. Install the ONNX support package:
   `Deep Learning Toolbox Converter for ONNX Model Format <https://www.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format>`_

   .. note::

      Support packages can be installed from MATLAB's HOME tab > Add-Ons >
      Get Add-Ons, then search for the package name.

3. Open MATLAB, navigate to the NNV root directory, and run:

   .. code-block:: matlab

      install

   .. note::

      If you restart MATLAB, rerun ``install.m`` or ``startup_nnv.m`` to
      restore paths. Alternatively, run ``savepath`` after installation
      (may require admin privileges).

Optional Packages
-----------------

**VGG16/VGG19 support** (for CNN verification on these architectures):

- `VGG16 <https://www.mathworks.com/matlabcentral/fileexchange/61733-deep-learning-toolbox-model-for-vgg-16-network>`_
- `VGG19 <https://www.mathworks.com/help/deeplearning/ref/vgg19.html>`_

**MATLAB Verification Library** (used in CAV 2023 comparison):

- `Deep Learning Toolbox Verification Library <https://www.mathworks.com/matlabcentral/fileexchange/118735-deep-learning-toolbox-verification-library>`_

**Framework converters** (load models from TensorFlow/PyTorch):

- `TensorFlow Converter <https://www.mathworks.com/matlabcentral/fileexchange/64649-deep-learning-toolbox-converter-for-tensorflow-models>`_
- `PyTorch Converter <https://www.mathworks.com/matlabcentral/fileexchange/111925-deep-learning-toolbox-converter-for-pytorch-models>`_

**Conformal Prediction** (requires Python 3.12+):

.. code-block:: bash

   # From the NNV root directory:
   python -m venv .venv
   # Windows: .venv\Scripts\activate
   # macOS/Linux: source .venv/bin/activate
   pip install -r requirement.txt

See :doc:`/user-guide/conformal-prediction` for detailed Python setup instructions.

.. _verify-installation:

Verify Installation
-------------------

Run the diagnostic tool to check your setup:

.. code-block:: matlab

   check_nnv_setup()

This checks:

- MATLAB version compatibility
- Required toolbox availability (critical, important, and optional tiers)
- Submodule paths (CORA, NNMT, HyST)
- Python environment for conformal prediction (if configured)

You should see a ``Ready!`` confirmation when all critical dependencies are satisfied.

You can also check the NNV version programmatically:

.. code-block:: matlab

   NNVVERSION()  % Returns 'NNV v3.0.0'

Execution Without Installation
------------------------------

NNV can be executed online without installing MATLAB through
`CodeOcean <https://www.codeocean.com>`_:

- **Latest** (CAV 2023): https://doi.org/10.24433/CO.0803700.v1
- CAV 2020 ImageStar: https://doi.org/10.24433/CO.3351375.v1
- CAV 2020 Tool paper: https://doi.org/10.24433/CO.0221760.v1

Uninstallation
--------------

To remove NNV from your MATLAB path:

.. code-block:: matlab

   cd code/nnv
   uninstall
