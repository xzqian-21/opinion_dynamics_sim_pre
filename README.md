# Opinion Polarization and Connected Disagreement — MATLAB Codes

This repository contains the MATLAB codes corresponding to the paper:

**Opinion polarization and its connected disagreement: Modeling and modulation**

The codes are used to generate all figures presented in the *main text* of the paper.

## Requirements
- MATLAB R2023b

## Usage Notes
1. Please add the `scripts/network_model` directory to the MATLAB path before running the codes.
2. All codes have been tested locally under MATLAB R2023b and are confirmed to run successfully.  
   The results from a representative single run are also included in this repository.
3. Some repeated numerical experiments may be computationally intensive and time-consuming on personal computers  
   (tested environment: MacBook Pro with Apple M1 chip and 8 CPU cores).  
   Users are advised to adjust the maximum number of workers in `parfor` loops according to their own hardware resources.

## Reproducibility
The scripts are designed to reproduce the numerical results and figures reported in the paper.
Users are encouraged to run the scripts in the order specified in the corresponding files.

## License
See the LICENSE file for details.

## Third-Party Code

This repository includes or adapts code from the following third-party sources:

- Portions of the code are adapted from external MATLAB implementations.
  The original license and copyright notices are retained in the corresponding source files. See scripts in `scripts/network_model`.
