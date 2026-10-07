# Changelog

All notable changes to the `twoDdipRY` package are documented in this file.

## [2026-10-07] - Output Beautification & Diagnostics Refactoring

### Added
- **Program Banner & Configuration Summary**:
  - Clear header banner identifying the `twoDdipRY` solver, physical model, and Rogers–Young (RY) closure.
  - Comprehensive **System Configuration & Parameters** block echoing grid dimensions ($N_r$, $r_{max}$, $q_{max}$), Hankel table status, restart files, interaction parameters ($\Gamma$, $z_1, z_2$, $\lambda_{12}$, $\sigma_{jk}$, $\gamma_{jk}$), solver tolerances, and state point scan ranges.
- **Structured State Point Reporting**:
  - Header banner for each evaluated $(\rho, x_2)$ state point.
  - **Thermodynamics & Equation of State**: Compressibility factor $Z = P/(\rho k_BT)$ with hard-disk ($Z_{HD}$) and dipolar ($Z_{dip}$) contributions, internal energy $U/(Nk_BT)$, and excess two-body entropy $S_2^{ex}/k_B$.
  - **Thermodynamic Consistency & Response**: Closure parameter $\eta$, inverse compressibility $\chi^{-1}$, virial pressure derivative $\partial P^*/\partial\rho$, isothermal compressibility $\chi_T = \partial\rho/\partial P^*$, and consistency discrepancy $f_{opt}$.
  - **Fluctuation & Spinodal Stability (Bhatia–Thornton)**: Zero-wavevector response matrix ($M_{rr}, M_{cc}, M_{rc}$), spinodal eigenvalues ($\lambda_1, \lambda_2$), physical stability classification (Stable, Near Spinodal Margin, Unstable Demixing), and concentration fluctuation $S_{cc}(0)$ with $1/S_{cc}(0)$.
  - **Structure Factor Highlights**: Peak of total structure factor $S(q)_{max}$, long-wavelength limit $S(0)$, fluctuation ratio $S(0)/S(q)_{max}$, disk size scaling ratio, and optical form-factor ratio $R_{01}/R_{02}$.
- **Multi-State Summary Table**:
  - Tabular summary of all scanned density and composition points upon execution completion.
  - Explicit list of generated output data files (`thermo.dat`, `gr.dat`, `sq.dat`, `sqopt.dat`, `srq.dat`, `solout.dat`).

### Changed
- **Picard Iteration Logging**:
  - Streamlined iteration logging for the physical state ($\rho$, `irho = 0`) with periodic progress updates every 50 iterations and a single summary line upon convergence.
  - Suppressed intermediate Picard iteration dumps for the finite-difference offset states ($\rho \pm \Delta\rho$), eliminating hundreds of lines of terminal noise.
- **Newton–Raphson Consistency Logging**:
  - Cleaned up output into single-line progress steps displaying $\eta$, $\chi^{-1}$, $\partial P^*/\partial\rho$, $f_{opt}$, and relative error.
  - Explicitly reports whether the closure was optimized via Newton–Raphson or evaluated with a fixed $\eta$.

### Fixed
- **Rogue `fort.3` Creation**:
  - Unit 3 (`2DdipdHD2c_out.dat`) was previously closed inside the inner loop, causing subsequent write operations on multi-state scans to write to an unformatted `fort.3` file. Unit 3 and unit 88 (`thermo.dat`) now remain open throughout the run and are properly closed at program termination.
- **Removed Debug Prints**:
  - Removed unformatted dumps of input parameters and raw density arrays (`rhot`, `rho`).
  - Removed leftover eigenvalue residual checks (`n1 = ...`, `n2 = ...`) and entropy sum diagnostics (`sums = ...`).
- **Safe Density Stepping**:
  - Added safety handling for `nrt = 0` to prevent division-by-zero during single-point density scans.
- **Repository Cleanliness**:
  - Added `results/` to `.gitignore` to prevent tracking large auxiliary binary and quadrature table files.
