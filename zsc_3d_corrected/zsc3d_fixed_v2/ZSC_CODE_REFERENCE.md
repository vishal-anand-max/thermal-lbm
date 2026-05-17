# ZSC Free-Energy LBM Code Reference Document
## For Vishal's Thermal LBM Solver Project

**Created:** February 28, 2026
**Purpose:** Comprehensive reference capturing detailed analysis of Arup Das's ZSC Fortran codes, identified bugs, fixes derived from Sudhakar's published work, and translation notes for JAX/XLB implementation. Upload this document at the start of any future conversation about ZSC multiphase implementation.

---

## 1. CODE INVENTORY

### Source: Arup Kumar Das, IIT Roorkee
Arup Das provided these codes to Vishal **before** Sudhakar (T. Sudhakar, NIT Uttarakhand) started his PhD. Sudhakar later used the same codebase (cleaned up) for his doctoral work on Taylor drop dynamics, publishing validated results in peer-reviewed journals.

### Files (4 Fortran codes + 2 input files):

| File | Description | Status |
|------|-------------|--------|
| `mainvelocity2profile_v2.f90` | 2D ZSC, dual-grid D2Q5/D2Q9 | Working, clean |
| `main_velocity_2profile.f90` | 2D ZSC, earlier version | Whitespace-only diff from v2 |
| `mainwalllatest_v2_forTecPlot_mod_April10.f90` | 3D ZSC, D3Q7/D3Q19, TecPlot output | Buggy BCs, many lines commented out |
| `main_walllatest.f90` | 3D ZSC, D3Q7/D3Q19, VTK output | More complete BCs but still has bugs |
| `properties.in` | Input parameters for 2D code | 1000:1 density ratio |
| `discrete.in` | Bubble positions for 2D code | 2 bubbles |

### Key relationship:
- `main_walllatest.f90` = more canonical 3D version (all BCs active, VTK output)
- `mainwalllatest_v2_forTecPlot_mod_April10.f90` = debugging branch (selectively disabled BCs, TecPlot output)
- Both 2D codes are functionally identical

---

## 2. ZSC MODEL PHYSICS (from ZSC 2006, JCP 218:353-371)

### 2.1 Two Distribution Functions

The ZSC model uses two populations:
- **f** (order parameter): Tracks the density difference φ = (ρ_H - ρ_L)/2, using D2Q5 (2D) or D3Q7 (3D)
- **g** (momentum): Tracks the average density n = (ρ_H + ρ_L)/2, using D2Q9 (2D) or D3Q19 (3D)

**NOTATION WARNING:** Sudhakar's paper uses `n` for average density and `φ` for density difference. Arup's Fortran codes use variable names like `phi_f`, `rhon`, `phin` which can be confusing. The physics is identical.

### 2.2 Chemical Potential

```
μ_φ = 4A·φ·(φ² - φ*²) - κ·∇²φ
```

Where:
- `φ* = (ρ_H - ρ_L)/2` (equilibrium order parameter)
- In Arup's code: `alpha4 = 4*alpha`, `phistar2 = phistar*phistar`
- So: `muPhi = alpha4*phin*(phin*phin - phistar2) - kappa*lapPhi`

### 2.3 Parameter Relationships

From surface tension σ and interface width W:

```
α = 0.75·σ / (W·φ*⁴)          [In code: alpha = 0.75*sigma/(IntWidth*phistar**4)]
κ = 0.5·α·(W·φ*)²             [In code: kappa = (IntWidth*phistar)**2 * alpha/2]
```

Or equivalently (from Sudhakar's paper eqs 1-2):
```
W = 2·√(2κ/A) / (ρ_H - ρ_L)
σ = 4·√(2κA)/3 · (ρ_H - ρ_L)³
```

Where A = α (amplitude parameter) and k = κ (gradient energy penalty).

### 2.4 Pressure Tensor (Full Thermodynamic Pressure)

From 2D code CollisionG (line 944):
```
p = α·(3φ⁴ - 2φ*²·φ² - φ*⁴) - κ_G·(φ·∇²φ + 0.5·|∇φ|²) + c_s²·ρ
```

Where `κ_G = κ/4` in the 2D dual-grid code (rescaled), or `κ_G = κ` in single-grid 3D.

### 2.5 Interfacial Force

```
F_x = μ_φ · ∂φ/∂x
F_y = μ_φ · ∂φ/∂y
F_z = μ_φ · ∂φ/∂z
```

These enter the momentum equation through the forcing term in the collision operator.

### 2.6 Mobility

```
θ_M = q·(τ_φ·q - 0.5)·δ·Γ
```

Where Γ (mobility coefficient) = 40 following ZSC 2006 recommendation.
- In the 2D code: `Gamma = 40000` (scaled by 1000 due to dual-grid factor)
- In 3D code: Γ should be ≈ 40

### 2.7 Equilibrium Distributions

**For f (order parameter, D3Q7):**
```
f_0^eq = -3·Γ·μ_φ + φ        [rest: Af0 = -3*Gamma*muPhi]
f_1-6^eq = 0.5·Γ·μ_φ          [moving: Af1 = 0.5*Gamma*muPhi]
```
With velocity terms: `f_i^eq += (φ/(2q)) · e_i · u` where q = 1/(τ_φ + 0.5)

**For g (momentum, D3Q19):**
```
g_0^eq = w_0 · (n - 6·φ·μ_φ)                    [Ag0]
g_1-6^eq = w_1 · (3·φ·μ_φ + n) + velocity terms  [Ag1]
g_7-18^eq = w_2 · (3·φ·μ_φ + n) + velocity terms  [Ag2]
```

**D3Q19 weights:** w_0 = 1/3, w_1 = 1/18, w_2 = 1/36

### 2.8 Streaming for f (Lax-Wendroff Style — UNIQUE TO ZSC)

The order parameter uses a modified streaming step, NOT standard LBM streaming:
```
f(ie,j,1,nxt) = η·f(i,j,1,now) + η2·f(ie,j,1,now)
```
Where:
- `η = 1/(τ_φ + 0.5)`
- `η2 = 1 - η`

This is visible in the Stream subroutine of both 2D and 3D codes.

### 2.9 Stencils

**3D Laplacian (isotropic, 18-neighbor):**
```
∇²φ = (1/36)·[Σ(12 diagonal neighbors) + 2·Σ(6 face neighbors) - 12·φ_center]
```
In code: `inv36 * (diag_sum + 2.0*(face_sum - 12*phi_center))`

**3D Gradient (12-neighbor):**
```
∂φ/∂x = (1/12)·[2·(φ_E - φ_W) + Σ(relevant diagonal differences)]
```

**2D Laplacian (9-neighbor, in dual-grid fine mesh):**
```
∇²φ = (1/6)·[Σ(4 diagonal) + 4·Σ(4 cardinal) - 20·φ_center]
```

**2D Gradient (8-neighbor):**
```
∂φ/∂x = (1/12)·[4·(φ_E - φ_W) + (diagonal differences)]
```

---

## 3. THE 2D CODE — DETAILED ARCHITECTURE

### 3.1 Dual-Grid Implementation

The 2D code uses a **dual-grid** approach (from ZSC 2006 paper):
- **Coarse grid** (momentum): D2Q9, size [xmin:xmax, ymin:ymax]
- **Fine grid** (order parameter): D2Q5, size [xmin_f:xmax_f, ymin_f:ymax_f] where dimensions are 2× the coarse grid

The fine grid provides improved interface resolution. All differential terms (Laplacian, gradients) are computed on the fine grid and transferred to the coarse grid.

### 3.2 Time Stepping

Per main iteration:
1. `CollisionG` — momentum collision+stream on coarse grid (1 full step)
2. `Update` — bilinear interpolation of u,p from coarse→fine grid
3. `CollisionF` → `Stream` → `Differentials` — order parameter half-step 1 on fine grid
4. `CollisionF` → `Stream` → `Differentials` — order parameter half-step 2 on fine grid
5. Swap `now`/`nxt` indices for coarse grid

Two fine-grid steps per one coarse-grid step maintains synchronization.

### 3.3 Update Subroutine (Bilinear Interpolation)

Maps coarse grid (i,j) to fine grid targets. For each coarse cell, 9 fine grid targets are updated:
- Direct copy at coincident nodes (targets 1,3,7,9)
- Average of 2 neighbors at edge midpoints (targets 2,4,6,8)
- Average of 4 neighbors at cell center (target 5)

### 3.4 Modifications from Clean ZSC (Research Hacks)

- **Line 483:** Parabolic initial y-velocity: `u_f(i,j,2) = -0.0D-5 * parabolic` (channel flow setup, coefficient essentially zero)
- **Lines 522-526:** Zou-He BCs for f at ymin/ymax walls
- **Lines 583-591:** Zou-He BCs for g at ymin/ymax walls
- **Line 927-931:** Conditional buoyancy: `if (phin > 0) sFy += 2·φ*·3.8E-5` (gravity on bubble phase only)
- **Line 959:** `delt` parameter (0.5 from properties.in) multiplies invTauRho in collision — sub-relaxation for stability
- **Commented-out bounce-back code** at domain edges — evidence of BC experimentation

### 3.5 What the 2D Code Models

**Two bubbles rising in a vertical channel with walls:**
- Domain: 160×320 (tall narrow channel)
- Two bubbles: (50,50,R=20) and (115,50,R=12)
- 1000:1 density ratio
- Zou-He walls at ymin/ymax, periodic in x
- Buoyancy-driven rise with bubble interaction

### 3.6 Validation Built In

- Laplace law: calculates Pin-Pout and compares to σ/R
- Mass conservation tracking
- Effective radius monitoring
- Convergence test on pressure difference (10 × tStats window)
- VTK output for Paraview visualization

---

## 4. THE 3D CODE — DETAILED ARCHITECTURE

### 4.1 Single-Grid Implementation

Unlike the 2D code, the 3D version uses a **single grid**:
- D3Q7 for order parameter (f), indices 0-6
- D3Q19 for momentum (g), indices 0-18
- Both stored as 5D arrays: `f(x,y,z,0:fdim,0:1)`, `g(x,y,z,0:gdim,0:1)`
- Double-buffered with `now`/`nxt` swapping

### 4.2 Time Stepping

Per main iteration:
1. `Collision` — handles BOTH f and g collision, Laplacian/gradient computation, macroscopic updates, AND wall BCs
2. `Stream` — f streaming (Lax-Wendroff) + wall BCs for f
3. Swap `now`/`nxt`

Note: The 3D code combines collision and macroscopic update in one monolithic subroutine, unlike the cleaner separation in 2D.

### 4.3 Geometry

The code attempts a **stepped channel** geometry:
- Outer walls: xmin, xmax, ymin, ymax (rectangular duct)
- Internal step walls at xmin+50, xmax-50, up to height zmax-50
- Creates narrow channel sections that open into a wider plenum
- Gravity with inclination: `gravX = f(cosθ)`, `gravZ = f(sinθ)`
- Eötvös number controls buoyancy/surface-tension ratio

### 4.4 What the 3D Code Models

**Bubble/drop rising through a stepped/L-shaped channel with tilted gravity** — precursor to Sudhakar's PhD work on Taylor drops in channel junctions (U-bends, T-junctions, sudden contractions).

---

## 5. BUGS IDENTIFIED IN THE 3D CODE

### Bug 1: Sequential Overwrite Wall BCs (CRITICAL)

**Location:** Init subroutine (lines 497-629), Collision subroutine, Stream subroutine — repeated everywhere walls are applied.

**The pattern:**
```fortran
g(i,j,k,1,now) = g(i,j,k,2,now)    ! Step A: dir 1 gets value of dir 2
g(i,j,k,2,now) = g(i,j,k,1,now)    ! Step B: dir 2 gets value of dir 1
                                      ! BUT dir 1 was already overwritten!
                                      ! Result: both directions = original dir 2 value
```

**What should happen (proper bounce-back):**
```fortran
temp = g(i,j,k,1,now)
g(i,j,k,1,now) = g(i,j,k,2,now)
g(i,j,k,2,now) = temp
```
Or equivalently, use a temporary array or apply bounce-back during streaming (the standard approach).

**Impact:** All wall boundary conditions destroy information. Both incoming and outgoing distributions end up with the same value, which is NOT bounce-back — it's closer to a symmetry/mirror condition. This may explain why the code partially "works" (symmetry BCs don't blow up) but gives incorrect physics at walls.

### Bug 2: Partially Disabled Inlet/Outlet BCs

**Location:** Init and Stream subroutines, the `c1`/`c2` buffer array sections.

In `mainwalllatest_v2_forTecPlot_mod_April10.f90`, most of the periodic-mapping lines at zmin (connecting narrow channel exits) are commented out with `!`. In `main_walllatest.f90`, they're active. Neither version appears to implement the **correct** scheme from Sudhakar's paper (Zou-He inlet + extrapolation+ghost outlet).

### Bug 3: Hardcoded Magic Numbers

Wall offsets (±50 from domain edges) are hardcoded throughout. The geometry should be parameterized.

### Bug 4: No Proper Contact Angle Implementation

Sudhakar's paper specifies 90° contact angle maintained by choosing k and A so μ_φ ≈ 0.25 at walls. The Fortran codes don't implement this — they just apply bounce-back (broken bounce-back at that) without wettability control.

### Bug 5: Periodic BC Code Commented Out

In Init (lines 434-440), the periodic boundary conditions are all commented out:
```fortran
!ni(xmin,:,:,2) = xmax
!ni(xmax,:,:,1) = xmin
...
```
The neighbor array is never corrected for boundary nodes, meaning streaming at boundaries accesses out-of-range indices (mitigated by the -1:+1 allocation padding, but physically wrong).

---

## 6. FIXES REQUIRED (Based on Sudhakar's Published Methodology)

### Published BC Scheme (from Sudhakar & Das, I&EC Research 2020):

1. **Walls:** Standard bounce-back with 90° contact angle (neutral wettability)
2. **Inlet:** Zou-He scheme with fully developed parabolic velocity profile for f_i's
3. **Outlet:** Constant pressure outlet using extrapolation + ghost layer
4. **Contact angle:** k and A chosen so μ_φ ≈ 0.25 at wall nodes

### Fix 1: Replace Wall BCs with Proper Bounce-Back

**For g (D3Q19):** During streaming, distributions that would stream into a wall node are reflected to the opposite direction:
```fortran
! At a wall node, AFTER streaming:
g(wall,opposite_dir,nxt) = g(wall,dir,nxt)  ! but must be done via pre/post-streaming logic
```

The cleanest implementation: apply bounce-back DURING the streaming step, not as a post-processing patch. For each wall node, the outgoing distributions become the incoming distributions in opposite directions.

**For f (D3Q7):** Same principle but with the Lax-Wendroff streaming. Need to ensure bounce-back is applied after the `η·f(i) + η2·f(neighbor)` combination.

### Fix 2: Implement Zou-He Inlet

At zmin (inlet plane), impose parabolic velocity profile:
```
u_z(x,y) = u_max · [1 - (2x/L_x - 1)²] · [1 - (2y/L_y - 1)²]
```

Apply Zou-He equations to determine unknown f_i's and g_i's from known velocity.

### Fix 3: Implement Pressure Outlet with Ghost Layer

At zmax (outlet), maintain constant pressure using:
- Extrapolate distributions from interior: `g(outlet) = 2·g(outlet-1) - g(outlet-2)`
- Or use ghost layer approach: add one extra layer beyond zmax with equilibrium distributions at target pressure

### Fix 4: Parameterize Geometry

Replace all hardcoded `±50` offsets with variables read from input file.

### Fix 5: Proper Neighbor Array at Boundaries

For wall boundaries, set neighbor array so that `ni(wall_node, into_wall) = wall_node` (self-reference), then handle bounce-back in collision/streaming.

---

## 7. SUDHAKAR'S PUBLISHED RESULTS (Validation Targets)

### Paper: Sudhakar & Das, Ind. Eng. Chem. Res. 2020, 59, 19045-19061

**System:** Kerosene Taylor drop in water medium through rectangular U-bend channel

**Parameters:**
- ρ_H = 1000 kg/m³ (water), ρ_L = 787 kg/m³ (kerosene)
- Channel cross-section: 1×1 cm² (a = 1 cm)
- Lattice resolution: 100×100 per cross-section
- Morton number: Mo = 4.2×10⁻¹⁰
- Eötvös number: Eo = 10
- Mach number: Ma < 0.01
- Time step: δt = 10⁻⁶ s (order of δ²)
- Mobility coefficient: Γ = 40
- Reynolds numbers tested: Re = 100, 1000, 1500

**Validation:** Rising kerosene drop in vertical rectangular channel compared against in-house experimental data (Figure 1f in paper). Good match in interface contour prediction.

**Key results for validation:**
- Taylor drop passage through U-bend takes ~5.5×10⁵ lattice time units at Re=1500
- Taylor drop passage takes ~7.4×10⁵ lattice time units at Re=1000
- Drop tip Reynolds number profiles (Figure 8)
- Averaged water film thickness at cross-sectional planes (Figures 5, 6)
- Drop elongation angle measurements (Figure 10)

**Physical scenarios simulated:**
1. Standard U-bend (Fig 1a) at Re=1000, 1500
2. U-bend with contraction at 2nd bend (Fig 1b)
3. U-bend with contraction at 1st bend (Fig 1c)
4. Narrow U-bend — half bend radius (Fig 1d)
5. Z-bend — opposite directional bends (Fig 1e)
6. Gravity orientation effects: upflow, downflow, horizontal

---

## 8. TRANSLATION NOTES FOR JAX/XLB IMPLEMENTATION

### 8.1 Architecture Mapping

| Fortran (Arup) | JAX/XLB (Your Implementation) |
|----------------|-------------------------------|
| BGK collision | ULBM cascaded collision |
| f array (D3Q7) | Order parameter lattice in XLB |
| g array (D3Q19) | Momentum lattice in XLB |
| Explicit loops | Vectorized JAX operations |
| Neighbor arrays `ni` | XLB's built-in streaming |
| Manual bounce-back | XLB's `HalfwayBounceBack` |
| Zou-He manual | XLB's `ZouHe` BC class |
| Sequential time-stepping | Functional JAX transformations |

### 8.2 Key Translation Challenges

1. **The Lax-Wendroff streaming for f** — this is NOT standard LBM streaming. XLB's default streaming won't work for the order parameter. You'll need a custom streaming function: `f_new = η·f_local + η2·f_neighbor`

2. **Dual distribution coupling** — g depends on φ (from f) through the chemical potential and interfacial force. This coupling must be maintained across JAX's functional paradigm.

3. **BGK → ULBM replacement** — The ZSC interfacial physics (μ_φ, pressure tensor, force) is independent of the collision operator. You replace ONLY the relaxation mechanism in g, keeping all ZSC terms identical. The f collision remains BGK-like (it's a convection-diffusion equation, not Navier-Stokes).

4. **Chemical potential computation** — Requires Laplacian of φ, which needs a 3D stencil convolution. In JAX: use `jax.scipy.signal.correlate` or manual stencil with `jnp.roll`.

### 8.3 Phase 4 Implementation Order

1. Get single-grid 3D ZSC working with BGK first (validate against Sudhakar)
2. Replace g collision with ULBM cascaded
3. Add thermal coupling (Phase 5)
4. Validate Laplace law, bubble rise, contact angle

### 8.4 Input Parameter Template

Based on Sudhakar's validated parameters:
```
rhoH = 1000.0    # Heavy fluid density
rhoL = 787.0     # Light fluid density
sigma = [from Eo] # Surface tension
IntWidth = 5.0    # Interface width in lattice units (typical)
Gamma = 40.0      # Mobility coefficient
tauRho = [from Mo] # Momentum relaxation time
tauPhi = 0.7      # Order parameter relaxation time (typical)
Eo = 10.0         # Eötvös number
```

---

## 9. RELATIONSHIP MAP

```
Arup Das (IIT Roorkee)
  ├── SOURCE of ZSC Fortran codes (2D + 3D)
  ├── Gave codes to Vishal (before Sudhakar's PhD)
  ├── Supervised Sudhakar's PhD using same codebase
  └── Contact for 3D ZSC debugging expertise

Sudhakar Thippavathini (NIT Uttarakhand → IIT Roorkee PhD)
  ├── PhD: "Understanding Interfacial Evolution of Taylor Drop in Liquid Filled Pipe Junction using LBM"
  ├── Published validated 3D ZSC results (I&EC Research 2020, Chem. Eng. Sci. 2018)
  ├── Used: ZSC + D3Q19 + BGK + bounce-back walls + Zou-He inlet + ghost-layer outlet
  └── Γ = 40, Eo = 10, Mo = 4.2e-10, Re = 100-1500

ZSC 2006 Paper (Zheng, Shu, Chew — JCP 218:353-371)
  └── Foundation for all codes above
```

---

## 10. FILE USAGE GUIDE

When starting a new conversation about ZSC implementation:

1. Upload this document (`ZSC_CODE_REFERENCE.md`)
2. Upload the relevant Fortran code(s) you want to work on
3. Optionally upload Sudhakar's paper + supplementary if detailed equation reference needed
4. Tell Claude: "I'm working on ZSC multiphase implementation. Read the reference document for context on Arup Das's codes and the fixes needed."

This document, combined with Claude's memory of the project context, should restore ~90% of the current conversation's knowledge.

---

*Document version 1.0 — February 28, 2026*
