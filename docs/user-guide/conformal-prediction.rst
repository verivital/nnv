Conformal Prediction Setup
==========================

.. rst-class:: lead

   NNV's conformal prediction (CP) module provides probabilistic verification
   with formal coverage guarantees. It requires a Python environment with
   PyTorch for training surrogate models.

----

Overview
--------

The CP verification method works by:

1. **Sampling** the input perturbation space
2. **Training a surrogate model** (linear or ReLU) in Python to approximate the
   network's behavior in specific directions
3. **Computing prediction regions** with conformal inference guarantees
4. **Verifying robustness** using the surrogate's prediction sets

This approach is particularly useful for networks where full reachability
analysis is intractable, such as large semantic segmentation networks.

Python Setup
------------

**Requirements:** Python 3.12+

Windows
^^^^^^^

.. code-block:: bash

   cd <nnv-root-directory>
   python -m venv .venv
   .venv\Scripts\activate
   pip install -r requirement.txt

macOS / Linux
^^^^^^^^^^^^^

.. code-block:: bash

   cd <nnv-root-directory>
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirement.txt

The ``requirement.txt`` includes:

- ``torch`` -- PyTorch for neural network training
- ``numpy`` -- Numerical computing
- ``scipy`` -- Scientific computing (for .mat file I/O)

Verify Setup
^^^^^^^^^^^^

After setup, verify in MATLAB:

.. code-block:: matlab

   check_nnv_setup     % Shows Python environment status

   % Or test directly
   python_path = cp_env()  % Returns path to Python executable

The ``cp_env()`` function automatically locates the virtual environment
at the NNV root directory.

Usage
-----

Basic CP robustness verification:

.. code-block:: matlab

   % Load network and create input set
   net = matlab2nnv(dlnet);
   IS = ImageStar(image, lb_pert, ub_pert);

   % Configure CP verification
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   reachOptions.coverage = 0.99;      % Coverage probability (default: 0.99)
   reachOptions.confidence = 0.99;    % Confidence level (default: 0.99)
   reachOptions.train_mode = 'Linear'; % Surrogate type: 'Linear' or 'ReLU'
   reachOptions.train_device = 'cpu';  % 'cpu' or 'gpu'
   reachOptions.train_epochs = 100;    % Training iterations
   reachOptions.train_lr = 0.001;      % Learning rate

   % Run verification
   target_class = 3;
   num_classes = 10;
   result = verify_robustness_cp(net, IS, reachOptions, target_class, num_classes);

Surrogate Model Types
^^^^^^^^^^^^^^^^^^^^^

- **Linear** (``train_mode = 'Linear'``): Faster training, simpler model. Good
  for networks with approximately linear behavior in the perturbation region.
- **ReLU** (``train_mode = 'ReLU'``): More expressive surrogate with ReLU
  activations. Better for highly nonlinear networks.

GPU Acceleration
----------------

To use GPU for surrogate model training:

1. Install CUDA-enabled PyTorch:

   .. code-block:: bash

      pip install torch --index-url https://download.pytorch.org/whl/cu121

2. Set the device in MATLAB:

   .. code-block:: matlab

      reachOptions.train_device = 'gpu';

Internal Architecture
---------------------

The CP verification pipeline uses these Python scripts located in
``code/nnv/engine/nn/Prob_reach/``:

- ``Trainer_Linear.py`` -- Linear surrogate model training
- ``Trainer_ReLU.py`` -- ReLU surrogate model training
- ``Direction_trainer.py`` -- Direction computation for dimensionality reduction

MATLAB orchestrates the pipeline:

1. ``CP_specification()`` -- Computes required sample sizes (N_dir, N, Ns)
   based on coverage and confidence parameters
2. ``Prob_reach()`` -- Main driver that calls Python trainers and builds
   the probabilistic reachable set
3. ``verify_robustness_cp()`` -- End-to-end robustness verification

Troubleshooting
---------------

**Python not found:**

- Ensure the ``.venv`` directory is in the NNV root (same level as ``requirement.txt``)
- Check that the Python executable exists:
  - Windows: ``.venv\Scripts\python.exe``
  - Unix: ``.venv/bin/python``

**Module not found errors:**

.. code-block:: bash

   # Activate venv and reinstall
   .venv/bin/activate  # or .venv\Scripts\activate on Windows
   pip install -r requirement.txt

**MATLAB can't find Python:**

The ``cp_env()`` function uses ``nnvroot()`` to locate the virtual environment.
Ensure NNV is properly installed (``install.m`` has been run).

See Also
--------

- :doc:`/theory/probabilistic` -- Mathematical foundations of conformal prediction
- `CP verification example <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/cifar10>`_
