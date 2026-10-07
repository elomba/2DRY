# 2D Dipolar Hard Disk Integral Equation Solver (`twoDdipRY`)

## 1. Overview
This package contains the liquid-state integral equation solvers for two-dimensional (2D) binary mixtures ($N_{sp} = 2$) of hard disks interacting through repulsive parallel dipole–dipole ($1/r^3$) interactions under the **Rogers–Young (RY)** closure.

The codes solve the 2D Ornstein–Zernike (OZ) relation using Hankel ($J_0$) transforms on Lado's quadrature grid, enforcing thermodynamic consistency between virial and compressibility routes via a 1D Newton–Raphson solver for parameter $\eta$. The latest production code also evaluates the zero-wavevector Bhatia–Thornton number–concentration fluctuation response matrix ($M_{rr}, M_{cc}, M_{rc}$) and computes spinodal stability eigenvalues ($\lambda_1, \lambda_2$).

---

## 2. Package Contents

### Source Codes
- **`2DdipRYN_final.f90`** (Recommended Production Version):
  - Fully optimized with Intel MKL BLAS `dgemv` for Hankel transforms (~3.5x faster end-to-end).
  - Precomputes invariant potentials and reciprocals in the iterative core.
  - Self-contained double-precision SLATEC module (`bessel_slatec_mod`) for modified Bessel functions.
  - Complete in-source documentation.
  - Structured, formatted terminal output with clear headers, real-time convergence tracking, and comprehensive thermodynamic and spinodal stability summaries.
- **`2DdipRYN2.f90`** (Historical version, Sep 2021): Regular density grid scan + spinodal fluctuation analysis.
- **`2DdipRYNfx.f90`** (Historical version, Jun 2021): Outer loop over discrete density list from `2DdipHD2c_list.dat`.
- **`2DdipRYN.f90`** (Historical base version, Jan 2021): Single-density evaluation with composition scan.

### Build System
- **`Makefile`**: Multi-compiler makefile supporting Intel (`ifx`, `ifort`) and GNU (`gfortran`).

### Input Files (Annotated)
- **`2DdipHD2c_map.dat`**: Primary input file for `2DdipRYN_final.f90` and `2DdipRYN2.f90`. Contains line-by-line inline explanations.
- **`2DdipHD2c_list.dat`**: Input file for `2DdipRYNfx.f90` (discrete density list mode).

### Documentation
- **`README.md`**: This guide.
- **`Changelog.md`**: Detailed log of updates, output beautification, and bug fixes.
- **`COMPARISON.md`**: Detailed comparative scientific and technical analysis of the codes, optimization benchmarks, and validation results.

### Reference Output Files
- **`thermo.dat`**: Reference thermodynamic and spinodal stability table.
- **`gr.dat`**: Reference radial distribution functions $g_{jk}(r)$.
- **`sq.dat`**: Reference structure factors $S_{jk}(q)$ and $S_{cc}(q)$.
- **`sqopt.dat`**: Reference optical/scattering folded structure factors.
- **`solout.dat`**: Reference restart solution vector $s_{SR}(r)$.

---

## 3. Compilation

### Using Intel oneAPI (Recommended)
Ensure the Intel compiler and MKL are in your environment (e.g. `module load intel` or `source /opt/intel/oneapi/setvars.sh`):

```bash
make
# or directly:
ifx -O3 -qmkl 2DdipRYN_final.f90 -o 2DdipRYN_final
```

### Using GNU Fortran (`gfortran`)
```bash
make FC=gfortran
# or directly:
gfortran -O3 2DdipRYN_final.f90 -lblas -llapack -o 2DdipRYN_final
```

---

## 4. Running the Code & Generating `ftable.dat`

The Hankel quadrature table file `ftable.dat` (~76 MB uncompressed) is **not included** in this archive to keep the package lightweight.

### First Run (Generating `ftable.dat`):
1. In `2DdipHD2c_map.dat`, verify line 1 has `newW = 1`:
   ```
   2500 1 0 1 0.00001 0.5 0.000001 0.00001 / ...
   ```
   The program will compute the Bessel roots and quadrature weights, and automatically save `ftable.dat`.

### Subsequent Runs:
2. Change `newW = 0` on line 1 of `2DdipHD2c_map.dat`:
   ```
   2500 1 0 0 0.00001 0.5 0.000001 0.00001 / ...
   ```
   The program will instantly load the precomputed `ftable.dat`.

### Execute:
```bash
./2DdipRYN_final
```

---

## 5. Input File Parameter Reference (`2DdipHD2c_map.dat`)

Each line corresponds to specific program variables:

| Line | Parameters | Description |
|:---:|---|---|
| **1** | `Nr iStart newW Gamma blend0 rmsMax rmscut` | Grid points ($N_r=2500$), restart flag (0=fresh, 1=restart), table flag (0=read, 1=generate `ftable.dat`), dipole coupling $\Gamma$, Picard blend, convergence thresholds |
| **2** | `nrt` | Number of density intervals for equispaced density scan ($inr = 0, \dots, nrt$) |
| **3** | `rtmin rtmax` | Total density scan range: $[\rho_{min}, \rho_{max}]$ |
| **4** | `nxf` | Number of composition points ($x_2$) |
| **5** | `xf(1:nxf)` | Mole fraction(s) of species 2 ($x_2$; species 1 has $1 - x_2$) |
| **6** | `z(1) z(2)` | Dipole moments / effective charges for species 1 and 2 |
| **7** | `Ncore(1,1) Ncore(1,2)` | Hard core grid radius indices for 1-1 and 1-2 pairs ($r \le r(N_{core})$ is inside core) |
| **8** | `Ncore(2,2)` | Hard core grid radius index for 2-2 pair |
| **9** | `'ftable.dat'` | Filename for Hankel transform tables |
| **10** | `'solin.dat'` | Input restart file (read only if `iStart = 1`) |
| **11** | `'solout.dat'` | Output restart file for converged solutions |
| **12** | `lambda shrink r01 R02` | 1-2 cross interaction scale ($\lambda_{12}$), core shrink factor, and disk scattering form-factor radii |
| **13** | `eta dsig tol osig` | Rogers–Young parameter $\eta$, Newton step $d\eta$, consistency tolerance, consistency flag (`.true.` = adjust $\eta$, `.false.` = fixed $\eta$) |

---

## 6. Output Files

- **`thermo.dat`** (Unit 88): Comprehensive thermodynamic summary table with a complete `#`-commented metadata header explaining every column, followed by 1-to-1 character-aligned columns (14 characters per column, right-aligned formatted as `(1x, 15(1x, f14.6))`):
  1. `rho_tot` : Total number density $\rho = \rho_1 + \rho_2$
  2. `x2` : Mole fraction of species 2 ($x_2 = \rho_2 / \rho_{tot}$)
  3. `U/NkT` : Reduced internal energy per particle $U / (N k_B T)$
  4. `P/(rho*kT)` : Virial compressibility factor $Z = P / (\rho k_B T)$
  5. `chi^-1` : Inverse isothermal compressibility from compressibility route ($\chi^{-1} = 1 - \hat{c}(0)$)
  6. `dP*/drho` : Virial pressure derivative with respect to density $\partial P^* / \partial \rho$
  7. `eta` : Rogers–Young closure consistency parameter $\eta$
  8. `Sq_max` : Maximum peak of total structure factor $S(q)_{max}$
  9. `Sq_0` : Zero-wavevector limit of total structure factor $S(q=0)$
  10. `Sq0/Sqmax` : Long-wavelength fluctuation ratio $S(0) / S(q)_{max}$
  11. `scale_21` : Effective core scaling factor $(S_{22}(0) / S_{11}(0))^{1/4}$
  12. `S2_ex/kB` : Excess two-body entropy per particle $S_2^{ex} / k_B$
  13. `lambda1` : Minimum eigenvalue of Bhatia–Thornton fluctuation matrix ($\lambda_1 \to 0$ signals spinodal boundary)
  14. `lambda2` : Maximum eigenvalue of Bhatia–Thornton fluctuation matrix
  15. `1/Scc(0)` : Inverse concentration fluctuation at origin (vanishes at critical demixing point)
- **`gr.dat`** (Unit 22): Radial distribution functions: column 1 is $r$, columns 2–5 are $g_{11}(r), g_{12}(r), g_{21}(r), g_{22}(r)$.
- **`sq.dat`** (Unit 23): Structure factors: wavevector $q$, reduced wavevector $q/\sqrt{\rho}$, $S_{jk}(q)$, total $S(q)$, form-factor folded structure factors, and concentration structure factor $S_{cc}(q)$.
- **`sqopt.dat`** (Unit 230): Form-factor folded optical structure factors.
- **`srq.dat`** (Unit 95): Scaled partial structure factors $(S_{jk}(q)\rho/\rho_j)$.
- **`solout.dat`**: Restart vector $s_{SR}(r; j, k)$ for subsequent runs with `iStart = 1`.

---

## 7. Terminal Output & Diagnostics

When executed, `2DdipRYN_final` prints a structured, ANSI syntax-colored report to standard output (all data files remain clean and free of escape codes):

1. **System Configuration**: Displays grid parameters ($N_r$, $r_{max}$, $q_{max}$), interaction parameters ($\Gamma$, $z_i$, $\lambda_{12}$, $\sigma_{jk}$, $\gamma_{jk}$), solver tolerances, and density/composition scan ranges highlighted in cyan and yellow.
2. **Convergence Progress**:
   - Displays periodic Picard iteration convergence (every 50 iterations) and final iteration count on the central physical state $\rho$ in cyan/yellow/white.
   - Displays step-by-step Newton–Raphson consistency progress ($\eta$, $\chi^{-1}$, $\partial P^*/\partial\rho$, consistency discrepancy $f_{opt}$, and relative error) with color-coded tags.
3. **State Point Results**:
   - **Thermodynamics & Equation of State**: Total compressibility factor $Z = P/(\rho k_BT)$ decomposed into hard disk ($Z_{HD}$) and dipolar ($Z_{dip}$) virials, reduced internal energy $U/(Nk_BT)$, and excess two-body entropy $S_2^{ex}/k_B$.
   - **Thermodynamic Consistency & Response**: Rogers–Young parameter $\eta$, inverse compressibility $\chi^{-1}$, virial pressure derivative $\partial P^*/\partial\rho$, isothermal compressibility $\chi_T = \partial\rho/\partial P^*$, and consistency discrepancy.
   - **Fluctuation & Spinodal Stability (Bhatia–Thornton)**: Response matrix ($M_{rr}, M_{cc}, M_{rc}$), spinodal eigenvalues ($\lambda_1, \lambda_2$), physical stability classification (`STABLE` in bold green, `NEAR SPINODAL MARGIN` in bold yellow, or `UNSTABLE / DEMIXING` in bold red), and concentration fluctuation $S_{cc}(0)$ and $1/S_{cc}(0)$.
   - **Structure Factor Highlights**: Total peak $S(q)_{max}$, zero-wavevector limit $S(0)$, ratio $S(0)/S(q)_{max}$, disk scaling ratio, and optical form-factor ratio $R_{01}/R_{02}$.
4. **Summary Table**: Colorized multi-column table summarizing all computed $(\rho, x_2)$ state points with their key thermodynamic and stability indices at program termination.

---

## 8. Numerical Monitoring & Fail-Safe Protection

To prevent runaway divergences from silently generating corrupted output files or ruining existing restart states, `2DdipRYN_final` incorporates automated monitoring via `Module monitor_mod`:

- **Real-Time IEEE Checks**: Continuously evaluates `ieee_is_finite` on all critical numerical quantities:
  - Picard residuals (`rms`), intermediate correlation arrays ($s_{SR}(r; j,k)$), and iteration limits (`iterMax`).
  - Virial pressures ($P_1, P_2$), internal energy ($U$), and compressibility integrals ($X$).
  - Bhatia–Thornton fluctuation matrix elements ($M_{rr}, M_{cc}, M_{rc}$), spinodal eigenvalues ($\lambda_1, \lambda_2$), and concentration fluctuation ($S_{cc}(0)$).
  - Newton–Raphson consistency derivatives ($\partial P^*/\partial\rho$, $f_{opt}$, $f'$), and parameter updates ($\eta, \tilde{\eta}$).
- **Pre-Flight Restart Protection**: Rigorously verifies that the indirect correlation array $s_{SR}$ and all state observables are strictly finite **before** opening or writing to `solout.dat`. If a NaN or Inf is detected, calculations halt immediately and previous valid restart files remain untouched.
- **Graceful File Cleanup**: In the event of a fatal numerical divergence or upon normal completion, `safe_close_all()` flushes write buffers and closes all active file handles (units `2`, `3`, `15`, `16`, `17`, `22`, `23`, `88`, `95`, `230`).
