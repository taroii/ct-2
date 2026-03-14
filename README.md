# CT-2: 3D CT/Tomosynthesis Image Reconstruction

Python + Fortran codebase for 3D computed tomography (CT) and tomosynthesis image reconstruction. Supports analytic (FDK) and iterative (ART, EM, Total Variation minimization) reconstruction methods.

## Setup

### 1. Create conda environment

```bash
conda create -n ct2 python=3.12
conda activate ct2
conda install -c conda-forge gfortran
pip install -r requirements.txt
```

### 2. Fix PATH for Fortran runtime (required on Windows)

**PowerShell:**
```powershell
$env:PATH = "$env:CONDA_PREFIX\Library\bin;$env:PATH"
```

**Git Bash / MSYS:**
```bash
export PATH="$CONDA_PREFIX/Library/bin:$PATH"
```

### 3. Compile Fortran extensions

From the `emil_code/` directory:

```bash
f2py -c totalvar.f -m totalvar
f2py -c ray_proj3D.f -m ray_proj3D
f2py -c tomo3DfortranAux.f -m tomo3Dfort
f2py -c view_proj3D.f -m view_proj3D
f2py -c view_proj3D_tomo.f -m view_proj3D_tomo
f2py -c radon3D.f -m radon3D
```

All 6 should produce `.pyd` (Windows) or `.so` (Linux/Mac) files in the `emil_code/` directory.

## Running examples

From the `emil_code/examples/` directory:

```bash
python fdk_head.py
```

This builds a head phantom, generates cone-beam projections, and reconstructs with the FDK algorithm.

### Other examples

| Script | Description |
|---|---|
| `fdk_head.py` | FDK reconstruction of head phantom (analytic, non-iterative) |
| `jawPhantomCCB.py` | Jaw phantom, circular cone-beam + TV minimization |
| `gridRandom10Views_eps_0p01.py` | TV reconstruction from 10 random views (linogram) |
| `gridRandom10VRays_eps_0p01.py` | TV reconstruction from random rays |
| `sphericalRadon.py` | Spherical Radon transform inversion with TV |
| `arttest.py` | Basic ART reconstruction test |

The iterative examples (jaw, grid, sphericalRadon) run until you write `stop` in their corresponding traffic file (e.g., `trafficA.txt`, `trafficC.txt`).

## Project structure

```
emil_code/
├── tomo3D.py                 # Core module: sinogram3D, linogram3D, image3D, phantom3D classes
├── FDK.py                    # Feldkamp-Davis-Kress analytic reconstruction
├── configs3D.py              # Sinogram geometry configs (circular cone-beam, Giotto, etc.)
├── linoConfigs3D.py          # Linogram geometry configs
├── sphericalRadonConfigs3D.py# Spherical Radon transform configs
├── Utils.py                  # Utility functions (MGH file I/O)
├── Makefile                  # Build targets for Fortran compilation
├── *.f                       # Fortran computational kernels
├── phantoms/                 # Test phantoms (head, jaw, breast, Defrise)
├── examples/                 # Reconstruction example scripts
├── tests/                    # Unit tests
└── materials_data/           # Material attenuation data
```
