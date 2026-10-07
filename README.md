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

- **`thermo.dat`** (Unit 88): Contains 15 columns:
  1. $\rho_{total}$ (total density)
  2. $x_2$ (mole fraction)
  3. $U / Nk_BT$ (internal energy)
  4. $P / (\rho k_BT)$ (virial compressibility factor)
  5. $\chi^{-1}(0)$ (inverse isothermal compressibility)
  6. $\partial P^* / \partial \rho$ (virial pressure derivative)
  7. $\eta$ (Rogers–Young parameter)
  8. $S(q)_{max}$ (maximum of total structure factor)
  9. $S(0)$ (structure factor at origin)
  10. $S(0) / S(q)_{max}$
  11. Packing ratio / scaling factor
  12. $S_2^{ex}$ (excess two-body entropy)
  13. $\lambda_1$ (minimum eigenvalue of fluctuation matrix; $\lambda_1 \to 0$ signals spinodal boundary)
  14. $\lambda_2$ (maximum eigenvalue of fluctuation matrix)
  15. $1 / S_{cc}(0)$ (inverse concentration fluctuation at origin; vanishes at demixing critical point)
- **`gr.dat`** (Unit 22): Radial distribution functions: column 1 is $r$, columns 2–5 are $g_{11}(r), g_{12}(r), g_{21}(r), g_{22}(r)$.
- **`sq.dat`** (Unit 23): Structure factors: wavevector $q$, reduced wavevector $q/\sqrt{\rho}$, $S_{jk}(q)$, total $S(q)$, form-factor folded structure factors, and concentration structure factor $S_{cc}(q)$.
- **`solout.dat`**: Restart vector $s_{SR}(r; j, k)$ for subsequent runs with `iStart = 1`.
