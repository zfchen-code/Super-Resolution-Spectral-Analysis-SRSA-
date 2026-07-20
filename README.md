# SRSA MATLAB Code and Reproducibility Examples

This package provides a compact MATLAB implementation of **Super-Resolution Spectral Analysis (SRSA)** and four examples accompanying the manuscript. The examples examine robustness to the assumed model order, resolution of closely spaced frequencies, resolution under Gaussian white noise, and application to a superconducting-gravimeter record. Precomputed MATLAB workspaces and editable MATLAB figures are included so that the reported results can be inspected without repeating the most computationally demanding calculations.

## Directory organization

```text

|-- README.md
|-- src/
|   `-- DoSRSA.m
|-- examples/
|   |-- Demo1_KRobustness.m
|   |-- Demo2_FrequencyResolution.m
|   |-- Demo3_NoiseResolution.m
|   `-- Demo4_M1RealDataApplication.m
|-- mat_file/
|   |-- Demo1_KRobustness.mat
|   |-- Demo2_FrequencyResolution_Rules.mat
|   |-- Demo2_FrequencyResolution_Truth.mat
|   |-- Demo3_NoiseResolution.mat
|   |-- Demo4_RealDataM1.mat
|   `-- m1_2011.mat
`-- figures/
    |-- Demo1_KRobustness_Figure (2).fig
    |-- Demo2_FrequencyResolution_Figure (1)_Rules.fig
    |-- Demo2_FrequencyResolution_Figure (1)_Truth.fig
    |-- Demo3_NoiseResolution_Figure (3).fig
    |-- Demo3_NoiseResolution_Figure (4).fig
    `-- Demo4_M1RealDataApplication_Figure (1).fig
```

- `src/` contains the SRSA implementation.
- `examples/` contains the simulation and application scripts. Each script describes the experiment, figure panels, computational environment, and representative runtime in its header comments.
- `mat_file/` contains the real-data input and saved workspaces. For Demo 2, the `*_Rules.mat` and `*_Truth.mat` files contain the results obtained with the rule-based and ground-truth model orders, respectively.
- `figures/` contains editable MATLAB figure files for the principal results.

## Software requirements

The code was developed and tested with **MATLAB R2024b**. It uses:

- MATLAB;
- Signal Processing Toolbox; and
- Parallel Computing Toolbox (`parfor` in `DoSRSA.m`).

The calculations do not require a GPU. Several examples repeatedly construct or apply Hankel operators and perform singular-value decompositions over multiple values of `Q_L`; consequently, a multicore CPU and sufficient RAM are recommended.

## Quick start

Start MATLAB in the `FinalCode` directory, or change to that directory, and add the source, example, and data directories to the MATLAB path:

```matlab
rootDir = pwd;  % rootDir must be the FinalCode directory
addpath(fullfile(rootDir,'src'));
addpath(fullfile(rootDir,'examples'));
addpath(fullfile(rootDir,'mat_file'));
```

Run an example by entering its script name, for example:

```matlab
Demo1_KRobustness
```

The `mat_file` directory must be on the MATLAB path when running Demo 4 because the script loads `m1_2011.mat`.

To inspect an editable figure without recomputation:

```matlab
openfig(fullfile(rootDir,'figures', ...
    'Demo4_M1RealDataApplication_Figure (1).fig'));
```

To inspect a saved workspace:

```matlab
workspaceFile = fullfile(rootDir,'mat_file','Demo4_RealDataM1.mat');
S = load(workspaceFile);
whos('-file',workspaceFile)
```

## Core function: `DoSRSA.m`

### Syntax

```matlab
[Bins,Q_LVEC,IM] = DoSRSA(CObsY,DeltaT,bins,K);
[Bins,Q_LVEC,IM] = DoSRSA(...,'Q_L',Q_L);
[Bins,Q_LVEC,IM] = DoSRSA(...,'AlgorithmMode',mode);
```

### Inputs

- `CObsY`: real- or complex-valued uniformly sampled time series. Row and column vectors are accepted, and all values must be finite.
- `DeltaT`: positive sampling interval. The time unit of `DeltaT` determines the frequency unit of `bins`; for example, seconds give hertz and years give cycles per year.
- `bins`: candidate frequencies at which the SRSA response is evaluated.
- `K`: assumed signal-subspace dimension, denoted by `TPG` inside the function and by `K_input`, `K_guess`, or `K_used` in the examples. It must be a positive integer smaller than the number of samples. For a real-valued sum of sinusoids away from zero frequency and the Nyquist frequency, each sinusoid generally contributes a conjugate pair of complex exponentials; the corresponding ideal subspace dimension is therefore twice the number of tones.
- `'Q_L'`: optional vector of Hankel row counts. Valid values satisfy `1 <= Q_L <= N-K`. Invalid requested values are removed with a warning. If this argument is omitted, the function selects at most 50 approximately evenly spaced integer values over the valid range.
- `'AlgorithmMode'`: either `'svd'` (default) or `'rsvd'`.

### Outputs

- `Bins`: output candidate-frequency vector. For real-valued input it equals `bins`. For complex-valued input the implementation returns the bilateral frequency grid in reversed-sign order.
- `Q_LVEC`: three-column array `[row_index, Q_L, L]`, where `Q_L` is the number of Hankel rows and `L = N-Q_L+1` is the number of Hankel columns, or window length.
- `IM`: row-wise min-max-normalized SRSA coherence response. Its rows correspond to `Q_LVEC`, and its columns correspond to `Bins`.

### Computational modes

- `'svd'` explicitly constructs each Hankel matrix and computes a complete singular-value decomposition. It is deterministic and serves as the reference high-accuracy mode, but has the highest computational and memory costs.
- `'rsvd'` uses a matrix-free randomized approximation of the dominant singular subspace. It is more efficient for large problems but may introduce a small approximation error. A fixed geometry-dependent random seed is used internally for reproducibility. If `K >= min(Q_L,L)`, the function falls back to the complete-SVD calculation.

### Minimal example

```matlab
DeltaT = 1;
t = (0:999)*DeltaT;
y = cos(2*pi*0.10*t) + 0.2*randn(size(t));
bins = 0:0.001:0.2;

[f,QL,response] = DoSRSA(y,DeltaT,bins,2, ...
    'Q_L',round(linspace(100,800,20)), ...
    'AlgorithmMode','svd');

imagesc(f,QL(:,2),response);
axis xy;
xlabel('Frequency (Hz)');
ylabel('Q_L');
colorbar;
```

## Interpreting the SRSA outputs

The fundamental SRSA result is the full two-dimensional response `IM(Q_L,f_c)`, rather than a curve at one arbitrarily selected value of `Q_L`. The parameter `Q_L` and the corresponding window length `L = N-Q_L+1` jointly determine the Hankel-matrix geometry. A useful interval must provide an appropriate compromise between frequency aperture and shifted-sample support.

In the examples, the displayed products have the following roles:

- **Full 2-D SRS:** the primary result used to inspect the evolution of candidate spectral features over the complete sampled `Q_L` range.
- **Local 2-D SRS:** a display window used only to improve readability within a stable-response interval; it does not replace or override the full 2-D result.
- **Diagnostic cross-section:** a one-dimensional response at a specified `Q_L`, used to illustrate the spectral shape for that particular Hankel geometry.
- **Final 1-D SRS:** a compact spectrum obtained by averaging responses over a specified stable `Q_L` interval.

A candidate peak is regarded as stable when it persists along the `Q_L` direction, retains a stable frequency position, and can be distinguished from localized background fluctuations or unstable responses near the boundaries of the sampled range.

SRSA peak height is a normalized frequency-coherence response determined by the noise-subspace projection. It is not a physical amplitude, power spectral density, or spectral energy. If physical or relative amplitudes are required, they should be estimated after frequency identification by fitting the identified components to the time-domain data or by applying an appropriate signal-subspace parameter-estimation method.

## Examples

### Demo 1: robustness to the assumed model order

`Demo1_KRobustness.m` constructs a noisy complex-valued time series containing four complex exponentials at positive and negative frequencies. The script evaluates diagnostic SRSA cross-sections at `Q_L = 1` and `L = N = 1000` for

```matlab
K_guess = [2,7,10,50,100,200];
```

It also calculates the full 2-D SRSA response using the known model order `K_true = 4` and displays a normalized bilateral Fourier amplitude spectrum for reference. This example demonstrates bilateral analysis of a complex-valued signal and shows how the principal frequency locations behave when the assumed model order is underestimated or conservatively overestimated.

The noise realization is random by default. Uncomment the supplied `rng(20260025,'twister')` line before running the script if an exactly repeatable realization is required.

Expected products:

- Figure (1): singular-value spectrum;
- Figure (2): six diagnostic cross-sections, the Fourier reference spectrum, and the full 2-D SRSA response.

The package includes `Demo1_KRobustness.mat` and `Demo1_KRobustness_Figure (2).fig`.

Representative runtime on the workstation specified in the script: approximately 1 min.

### Demo 2: resolution of closely spaced frequencies

`Demo2_FrequencyResolution.m` analyzes two noise-free, real-valued multitone signals containing five and nine closely spaced sinusoids. Their ground-truth signal-subspace dimensions are `K1_true = 10` and `K2_true = 18`. Because the true model order is generally unknown in applications, the distributed calculation uses

```matlab
K1_input = 180;
K2_input = 180;
```

These practical values are selected from the elbow regions of the corresponding singular-value spectra. The ground-truth alternatives are documented in the script and are represented by the supplied `*_Truth` workspace and figure.

The example provides the simulated series, singular-value spectra, Fourier reference spectra, full 2-D SRSA responses, local 2-D display windows, and final stacked 1-D SRS results. The local displays use the common interval `Q_L in [8000,10000]`; the final 1-D spectra are constructed from the explicitly specified `stackRows` in the script. The complete 2-D responses remain the basis for interpreting the spectral results.

This is a large calculation based on daily sampling over approximately 41 years and dense candidate-frequency grids. Complete-SVD mode can require substantial memory and computation time.

Supplied products:

- `Demo2_FrequencyResolution_Rules.mat` and `Demo2_FrequencyResolution_Figure (1)_Rules.fig`;
- `Demo2_FrequencyResolution_Truth.mat` and `Demo2_FrequencyResolution_Figure (1)_Truth.fig`.

Representative runtime for the rule-based complete-SVD calculation on the workstation specified in the script: approximately 25 min.

### Demo 3: two-tone resolution in Gaussian white noise

`Demo3_NoiseResolution.m` evaluates two-tone resolution over a nominal 50-year interval. Each time series contains unit-amplitude components at `3 cpy` and `3 + Delta_f cpy`, with

```matlab
Delta_f = [0.040,0.020,0.0133];  % cpy
```

The phases are selected so that the two components are in phase at the center of the observation interval. Nine groups are formed by combining the three frequency separations with the following pointwise noise-level indices:

```matlab
SNR9999Grid = [1,  5,  9; ...
               5,  9, 13; ...
               9, 13, 17];
```

The Gaussian white-noise standard deviation is calibrated from the one-sided Fourier-amplitude normalization. At a fixed nonzero frequency bin, the pointwise 99.99th-percentile Rayleigh noise amplitude is set to `1/SNR9999` of the theoretical amplitude of a single tone.

For each of the nine groups, the script generates 100 independent Monte Carlo realizations using

```matlab
rng(20260025,'twister');
```

It then computes singular-value spectra, Fourier amplitude spectra, and SRSA responses for all 900 realizations. The SRSA calculations use matrix-free randomized SVD. Ensemble averaging is applied to reduce realization-specific noise fluctuations and to test whether the two candidate frequencies remain consistently distinguishable.

Expected products:

- Figure (1): singular-value spectra for the 100 realizations in each group;
- Figure (2): Fourier amplitude spectra for all realizations, together with the prescribed noise-amplitude threshold and theoretical single-tone amplitude;
- Figure (3): full ensemble-mean 2-D SRSA responses for all nine groups;
- Figure (4): local 2-D SRSA displays over `Q_L in [2500,3500]` and comparisons between the stacked 1-D SRS and ensemble-mean Fourier spectra.

This is the most computationally intensive example. It performs 900 singular-value analyses and 900 SRSA evaluations. The supplied `Demo3_NoiseResolution.mat` workspace is approximately 1.45 GB; sufficient storage and memory should be available before loading it.

The package includes `Demo3_NoiseResolution.mat`, `Demo3_NoiseResolution_Figure (3).fig`, and `Demo3_NoiseResolution_Figure (4).fig`.

Representative runtime for the complete experiment on the workstation specified in the script: approximately 70 min.

### Demo 4: application to the 2011 M1 superconducting-gravimeter record

`Demo4_M1RealDataApplication.m` demonstrates the data-driven SRSA workflow in the target frequency band of the Earth's fundamental spheroidal mode `0S2`. The input file `m1_2011.mat` contains two versions of the M1 record in `series_m1`: the original preprocessed and model-corrected record, sampled at 60 s, and a band-pass-filtered and downsampled series, sampled at 300 s.

The script first compares the Fourier power spectra of the two series to confirm that the principal spectral structure within the target band is retained after preprocessing. It then examines the singular-value spectrum of the 300-s series and evaluates SRSA for three adjacent model orders selected near its elbow:

```matlab
K_used = [50,55,60];
```

For each model order, the script displays the full frequency-`Q_L` response, the stable interval `Q_L in [641,2321]`, and the corresponding mean-stacked 1-D SRS. Published normal-mode frequencies are plotted only as interpretive references; they are not used as fixed-frequency constraints, peak-picking criteria, or inputs to the SRSA calculation.

The package includes `Demo4_RealDataM1.mat` and `Demo4_M1RealDataApplication_Figure (1).fig`.

Representative runtime on the workstation specified in the script: approximately 1 min.

## Reproducibility and computational notes

1. Run the scripts without modification to reproduce the stated numerical experiments. Changing the frequency-grid density, sampled `Q_L` values, model order, Monte Carlo count, random seed, or decomposition mode defines a different numerical experiment.
2. When `'Q_L'` is not supplied, `DoSRSA.m` evaluates at most 50 approximately evenly spaced valid Hankel geometries. Supplying `'Q_L'` provides explicit control and can reduce the cost of exploratory calculations.
3. Complete-SVD results are deterministic. The randomized-SVD implementation uses fixed internal seeds so that repeated calls with identical inputs and settings are reproducible within normal floating-point tolerances.
4. The scripts use `parfor` to distribute calculations over `Q_L`. Runtime depends on the number of workers, processor performance, memory bandwidth, candidate-frequency-grid size, and the number and dimensions of the Hankel operators.
5. Explicit Hankel matrices and complete singular-value decompositions can require substantial memory. For exploratory large-scale calculations, reduce the number of candidate frequencies or sampled `Q_L` values, or use `'rsvd'` when its approximation is acceptable. Any such change should be reported because it modifies the published experiment.
6. The distributed scripts create MATLAB figures in the current session but do not automatically overwrite the supplied `.mat` or `.fig` files.
7. MATLAB `.fig` files can be opened with `openfig`. Use MATLAB export functions when a standard raster or vector format is required.

## Citation and contact

If this code contributes to published work, please cite the associated SRSA manuscript and identify the code version used.

Author: Zhifeng Chen  
E-mail: zfchen@whu.edu.cn
