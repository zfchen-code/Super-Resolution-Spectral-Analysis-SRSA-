#  Super-Resolution Spectral Analysis (SRSA) – MATLAB Implementation

This repository provides a MATLAB implementation of **SRSA (Super-Resolution Spectral Analysis)** together with two demonstration scripts:

1. **High-frequency resolution demo** – illustrates the ability of SRSA to resolve closely spaced spectral components under noise-free conditions.
2. **Robustness demo** – demonstrates the insensitivity of SRSA to the initial model-order choice (`K_used`) and compares the results with the Fourier spectrum.

---

## Author
**Zhifeng Chen**  
Email: zfchen@whu.edu.cn  
Tested platform: MATLAB R2024b  
Created: 2026-01-20

---

## Repository Structure

### Core function
- `DoSRSA.m`  
  Core SRSA routine for computing the 2D super-resolution spectrum (interference-response matrix) by scanning window length \(L\) and candidate frequencies.

### Demonstration scripts
- **Demo 1 – High-frequency resolution (noise-free)**  
  Shows that SRSA can resolve closely spaced frequencies beyond the classical Fourier resolution limit.

- **Demo 2 – Robustness to K_used**  
  Shows that SRSA remains stable for a wide range of model-order guesses and compares results with FFT.

---

## Quick Start

### Requirements
- MATLAB R2024b (tested)
- No additional toolboxes required (only standard functions such as `svd`, `hann`, `hankel`)

### Setup
Make sure `DoSRSA.m` is on the MATLAB path.

### Run the demos
Simply execute the corresponding scripts in MATLAB:

```matlab
% High-frequency resolution test
run Demo1.m

% Robustness test
run Demo2.m
