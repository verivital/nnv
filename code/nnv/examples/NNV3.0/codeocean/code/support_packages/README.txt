MATLAB Support Packages for NNV 3.0
====================================

To enable ONNX-based tests (FairNNV, ProbVer, VideoStar), you need to create
and upload sppFile.zip containing the ONNX converter support package.

Instructions:
-------------

1. Open MATLAB (R2024a or compatible version)

2. Install the ONNX converter support package:
   - Go to Add-Ons > Get Add-Ons
   - Search for "Deep Learning Toolbox Converter for ONNX Model Format"
   - Click Install

3. Run the create_sppFile.m script (in the codeocean folder):
   >> cd('/path/to/codeocean')
   >> create_sppFile

4. Upload the generated sppFile.zip to this folder on CodeOcean:
   /code/support_packages/sppFile.zip

   OR to:
   /data/support_packages/sppFile.zip

5. Re-run the CodeOcean capsule

Note: The support package must match the MATLAB version used on CodeOcean.
CodeOcean uses R2024a, so create sppFile.zip using MATLAB R2024a if possible.
