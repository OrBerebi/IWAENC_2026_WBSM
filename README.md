# Field-of-View Informed Binaural Signal Matching for Head-Worn Arrays

This repository contains the evaluation and simulation codebase for the paper **"Field-of-View Informed Binaural Signal Matching for Head-Worn Arrays"** (originally developed under the working title *Weighted Binaural Signal Matching for Improved Binaural Reproduction Within a Constant Field of View*).

It provides the MATLAB simulation framework used to evaluate the FoVi-BSM method against the standard signal-independent Binaural Signal Matching (BSM) and the signal-dependent COMPASS-BSM baseline under complex acoustic conditions, specifically utilizing a wearable glasses microphone array (Project Aria prototype).

## Overview

The codebase is structured around the two primary numerical evaluations presented in the paper:
1. **Evaluation 1 (Anechoic Conditions):** Evaluates spatial (ILD, ITD) and spectral (NMSE, MagErr, BSD) reproduction errors for single far-field sources, including the parameter search for the spatial weighting constant $\beta$.
2. **Evaluation 2 (Reverberant Conditions):** A Monte Carlo simulation framework utilizing the Image Source Method to test multi-source acoustic scenarios across varying Direct-to-Reverberant Ratios (DRR) and room sizes.

## Repository Contents

The root directory contains the core execution scripts and dependency folders:

* **`startup_script.m`**: Initialization script to add all necessary folders and third-party tools to your MATLAB path.
* **`main.m`**: Main execution file.
* **`anechoic_analysis_FOV_BSM_v5.m`**: Script to reproduce Evaluation 1 (Anechoic baseline comparison and $\beta$ ablation study).
* **`main_monte_carlo.m`**: Script to reproduce Evaluation 2 (Monte Carlo simulations across reverberant room conditions).
* **`bsm_binaural_renderer_V3.m`**: Core function for rendering the binaural signals utilizing the BSM and FoVi-BSM filters.
* **`room_auralisation_v4.m`**: Handles the acoustic environment simulation and generation.

**Data & Function Folders:**
* `ATF/` and `HRTF/`: Contain the measured array transfer functions (Project Aria) and reference dummy head HRTFs.
* `dry_signals/`: Contains the anechoic source stimuli (e.g., TIMIT speech samples).
* `+BSM_toolbox/`, `+RoomParams/`, `+image_method/`: Object-oriented MATLAB packages for the core algorithms and room modeling.
* `aria_functions/`, `monte_carlo_functions/`, `parametric_encoding_functions/`, `room_auralization/`, `local_functions/`: Sub-modules containing the specific mathematical implementations and baseline metrics (including the parametric COMPASS-BSM rendering).
* `MUSHRA/`: Assets related to the subjective listening test.

## Prerequisites

-   **MATLAB** (R2022b or newer recommended)
-   **Signal Processing Toolbox**

## Usage

1.  **Initialize:** Open MATLAB and run `startup_script.m` to ensure all required directories (`+BSM_toolbox`, `parametric_encoding_functions`, etc.) are in your path.
2.  **Anechoic Evaluation:** Run `anechoic_analysis_FOV_BSM_v5.m` to generate the spectral and spatial cue error metrics over the dense evaluation grid.
3.  **Reverberant Evaluation:** Run `main_monte_carlo.m` to execute the multi-source reverberant evaluations across the simulated room conditions.

## Citation

If you use this code or our findings in your research, please cite our paper:

```bibtex
@article{berebi2026fovibsm,
  title={Field-of-View Informed Binaural Signal Matching for Head-Worn Arrays},
  author={Berebi, Or and Ben-Hur, Zamir and Alon, David Lou and Rafaely, Boaz},
  journal={TBD},
  year={2026}
}
