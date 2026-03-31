VNN-COMP & ARCH-COMP
====================

.. rst-class:: lead

   NNV participates in two major annual verification competitions,
   contributing to community benchmarking and tool evaluation.

----

VNN-COMP (Verification of Neural Networks Competition)
------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 10 50 40

   * - Year
     - Contribution
     - Source
   * - 2025
     - Star set reachability approach (competition contribution)
     - `SAIV 2025 <https://doi.org/10.1007/978-3-031-99991-8_15>`_
   * - 2024
     - Full benchmark participation with multi-solver fallback
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/VNN_COMP2024>`_
   * - 2023
     - `VNN-COMP 2023 <https://arxiv.org/abs/2312.16760>`_
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/VNN_COMP2023>`_
   * - 2021
     - `VNN-COMP 2021 <https://arxiv.org/abs/2109.00498>`_
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/VNN_COMP2021>`_

NNV will participate in both **VNN-COMP 2026** and **ARCH-COMP 2026**.

**Competition strategy**: NNV uses a three-phase approach:

1. Random falsification (100 samples) for quick counterexample search
2. Reachability analysis (approx-star) if no counterexample found
3. Multiple LP solver fallback (Gurobi → GLPK → linprog)

ARCH-COMP (Applied Verification of Continuous and Hybrid Systems)
-----------------------------------------------------------------

NNV participates in the **AINNCS** (Artificial Intelligence and Neural Network
Control Systems) category:

.. list-table::
   :header-rows: 1
   :widths: 10 50 40

   * - Year
     - Report
     - Source
   * - 2025
     - `ARCH-COMP25 AINNCS Report <https://mediatum.ub.tum.de/doc/1839200/nof8duj8633x0tjn7zvu9rm1v.ARCH25_pages_71-121_AINNCS.pdf>`_
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/ARCH_COMP2025>`_
   * - 2024
     - `ARCH-COMP24 AINNCS Report <https://easychair.org/publications/paper/WsgX>`_
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/ARCH_COMP2024>`_
   * - 2023
     - `ARCH-COMP23 AINNCS Report <https://easychair.org/publications/paper/Vfq4b>`_
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/ARCH_COMP2023>`_
   * - 2022
     - `ARCH-COMP22 AINNCS Report <https://easychair.org/publications/paper/C1J8>`_
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/ARCH_COMP2022>`_
   * - 2021
     - `ARCH-COMP21 AINNCS Report <https://easychair.org/publications/paper/Jq4h>`_
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/ARCH_COMP2021>`_
   * - 2020
     - `ARCH-COMP20 AINNCS Report <https://easychair.org/publications/paper/Jvwg>`_
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/ARCH_COMP2020>`_
   * - 2019
     - Initial participation
     - `Submission code <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/ARCH_COMP2019>`_

Other Conference Submissions
----------------------------

NNV has reproducible code for publications at:

.. list-table::
   :header-rows: 1
   :widths: 15 35 50

   * - Venue
     - Year(s)
     - Topics
   * - CAV
     - 2020, 2021, 2023
     - Core tool papers, ImageStar, semantic segmentation
   * - FM
     - 2019
     - Star-based reachability for DNNs
   * - FormaliSE
     - 2023, 2024, 2025
     - BNNs, malware detection, video verification
   * - FMICS
     - 2023
     - Time-dependent neural networks
   * - ICAIF
     - 2024
     - FairNNV fairness verification
   * - HSCC
     - 2022, 2023
     - RNN verification, evaluation methods
   * - NeurIPS
     - 2025
     - Probabilistic semantic segmentation verification

All submission code is available in
`examples/Submission/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission>`_.
