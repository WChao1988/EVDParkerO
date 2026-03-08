# Enhanced Varying-Density Parker-Oldenburg Inversion for Moho Depth and Seafloor Topography

This repository provides Python implementations of the **Parker–Oldenburg method** for estimating the depth of two density interfaces from gravity anomalies:

- **Moho depth** – class `EVDMOHO` in `faatomoho.py`
- **Seafloor topography** – class `EVDBATHYMETRY` in `faatobathymetry.py`

Both algorithms iteratively solve the nonlinear inverse problem in the frequency domain using a series expansion. The code is designed for reproducibility and can be used to generate the results presented in [The Enhanced Varying-Density Parker-Oldenburg
Method for Moho and Seafloor Imaging] (if applicable).

---

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
- [Key Parameters](#key-parameters)
- [Usage Examples](#usage-examples)
  - [Moho Inversion](#moho-inversion)
  - [Bathymetry Inversion](#bathymetry-inversion)
- [Reproducing Paper Results](#reproducing-paper-results)
- [Citation](#citation)
- [License](#license)
- [Zenodo Archive](#zenodo-archive)

---

## Requirements

- Python 3.6 or higher
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/) (specifically `scipy.fftpack` and `scipy.signal.windows`)
- [abc](https://docs.python.org/3/library/abc.html) (specifically `ABCMeta` and `abstractmethod`)
- [math](https://docs.python.org/zh-cn/3/library/math.html) (specifically `factorial`)


---

## Installation

Clone the repository and install the required packages:

```bash
git clone https://github.com/yourusername/parker-oldenburg-inversion.git
cd parker-oldenburg-inversion

---
```
## Data Preparation

Input data are expected as **2D NumPy arrays saved in text format** (e.g., `.txt` files readable by `numpy.loadtxt`).  
The repository includes two demo scripts that expect the following folder structure:

```bash
data_require_moho_inversion/
    gravity_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt
    water_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt
    sediment_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt
    crust_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt
    moho_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt

data_require_bathymetry_inversion/
    gravity_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt
    water_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt
    sediment_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt
    longwave_topo_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt
```


- **`gravity`** : free‑air gravity anomalies (mGal)
- **`water`** : water depth / bathymetry (km, positive downward)
- **`sediment`** : sediment thickness (km)
- **`crust`** : crustal thickness (km) – used only for Moho inversion
- **`moho`** : Moho depth (km) – used as reference / target in Moho inversion
- **`longwave_topo`** : long‑wavelength topography correction (km) – used only in bathymetry inversion

All matrices must have the **same dimensions** (`nrow`, `ncol`).  
The physical dimensions of the grid (`longrkm`, `longckm`) are the total lengths in the row and column directions (in **km**) – these are used to compute wavenumbers.

> **Note:** The filenames encode the geographic bounds (`lat_up`, `lat_down`, `lon_left`, `lon_right`). Adjust these placeholders according to your actual data coverage.

##Key Parameters
| Parameter            | Description                                                                                              | Typical Value / Default          |
|----------------------|----------------------------------------------------------------------------------------------------------|----------------------------------|
| `grid_size`          | Dimensions of the input gravity grid (rows, cols) – automatically inferred from `faa_matrix.shape`.     | –                                |
| `padding`            | Not explicitly used; edge effects are mitigated by a **Tukey window** (`edge=0.02` in `twkey` method).  | Tukey window `edge=0.02`         |
| `fft_handling`       | FFT performed with `scipy.fftpack.fft2` / `ifft2`. Data are multiplied by a Tukey window before FFT.    | –                                |
| `stopping_criteria`  | Currently the iteration runs for a fixed number of steps (`t`). Convergence can be implemented manually.| –                                |
| `max_iterations`     | Number of iterations (`t` parameter in `downward` method).                                              | 4–8 (examples)                   |
| `wh`                 | Cutoff frequency (rad/km) for the low‑pass filter.                                                      | 0.2 (examples)                   |
| `alpha`              | Exponent controlling the filter roll‑off (high values = sharper cutoff).                                 | 8–64 (examples)                  |
| `delta_rhos`         | Density contrasts (g/cm³) for each interface. Dictionary with keys 0 and 1.                             | Moho: `{0:2.5, 1:1.63}`<br>Bathymetry: `{0:1.3, 1:0.85}` |
| `reference_depths`   | Reference depths (km) for each interface. Dictionary with keys 0 and 1.                                 | Moho: `{0: -moho_mean, 1: -crust_mean}`<br>Bathymetry: `{0: -(water-sediment).mean(), 1: -water.mean()}` |
| `mus`                | Decay factors for the series expansion. Typically set to small values (0–0.01).                         | `{0:0, 1:0.01}` (examples)       |
| `delta_rho_initial`  | Initial density contrast for bathymetry inversion (g/cm³).                                              | 1.77 (example)                   |
| `longwave_topo`      | Long‑wavelength topography correction (km) – used only in bathymetry inversion.                         | Provided as a matrix              |
| `delta_bnds`         | Initial interface undulations (km) for Moho inversion. Dictionary with keys 0 and 1.                    | `{0: None, 1: (water-sediment).mean()-removed}` |

> **Note:** The `downward` method runs for a fixed number of iterations (`t`). If you need to implement a residual‑based stop, you can modify the loop to break when the change between iterations falls below a threshold (e.g., 1e-6).

## Usage Examples

### Moho Inversion
```python
import numpy as np
from faatomoho import EVDMOHO

# Geographic bounds of the study area
lat_up, lat_down, lon_left, lon_right = 20, 15, 129, 134

# Load data
gravity = np.loadtxt(f'./data_require_moho_inversion/gravity_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
water   = np.loadtxt(f'./data_require_moho_inversion/water_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
sediment= np.loadtxt(f'./data_require_moho_inversion/sediment_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
crust   = np.loadtxt(f'./data_require_moho_inversion/crust_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
moho    = np.loadtxt(f'./data_require_moho_inversion/moho_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')

# Grid dimensions (km)
longrkm = 555.9746332227937
longckm = 530.2272075400947

# Parameter dictionary
params = {
    'faa_matrix': gravity - gravity.mean(),          # remove mean
    'delta_bnds': {0: None,
                   1: (water - sediment) - (water - sediment).mean()},
    'delta_rhos': {0: 2.5, 1: 1.63},
    'mus': {0: 0, 1: 0.01},
    'reference_depths': {0: -moho.mean(), 1: -crust.mean()},
    'longrkm': longrkm,
    'longckm': longckm,
    'wh': 0.2,
    'alpha': 64
}

# Run inversion
solver = EVDMOHO(**params)
moho_undulation = solver.downward(t=4)                # 4 iterations
moho_depth = -params['reference_depths'][0] + moho_undulation   # positive depth

np.save('moho_depth.npy', moho_depth)
```
### Bathymetry Inversion

```python

import numpy as np
from faatobathymetry import EVDBATHYMETRY

# Geographic bounds
lat_up, lat_down, lon_left, lon_right = 19, 16, 130, 133

# Load data
gravity = np.loadtxt(f'./data_require_bathymetry_inversion/gravity_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
water   = np.loadtxt(f'./data_require_bathymetry_inversion/water_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
sediment= np.loadtxt(f'./data_require_bathymetry_inversion/sediment_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
longwave_topo = np.loadtxt(f'./data_require_bathymetry_inversion/longwave_topo_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')

# Grid dimensions (km)
longrkm = 333.5847799336763
longckm = 318.1421721146795

# Parameter dictionary
params = {
    'faa_matrix': gravity,
    'sediment_matrix': sediment,
    'delta_rhos': {0: 1.3, 1: 0.85},
    'delta_rho_initial': 1.77,
    'mus': {0: 0.01, 1: 0.01},
    'reference_depths': {0: -(water - sediment).mean(), 1: -water.mean()},
    'longrkm': longrkm,
    'longckm': longckm,
    'wh': 0.2,
    'alpha': 8,
    'longwave_topo': longwave_topo
}

# Run inversion
solver = EVDBATHYMETRY(**params)
bathymetry_undulation = solver.downward(t=6)          # 6 iterations
bathymetry_depth = -params['reference_depths'][1] + bathymetry_undulation   # positive depth

np.save('bathymetry_depth.npy', bathymetry_depth)
```

##Reproducing Paper Results
To reproduce the figures and numerical results presented in the associated paper:

1.Download the input datasets from [DOI/data-link] (if not already provided in the repository).

2.Place the data files in the appropriate folders: data_require_moho_inversion/ and data_require_bathymetry_inversion/.

3.Run the demo scripts:
```bash
python demo_moho.py
python demo_bathymetry.py
```
The output depth grids `(.npy files)` will be saved in the current directory.
If you have modified any parameters, please refer to the Key Parameters section for default values used in the paper.

##License
This project is licensed under the MIT License – see the LICENSE file in the repository for full details.
##Zenodo Archive
A permanent snapshot of this repository is archived on Zenodo, with the DOI badge linked below:
https://zenodo.org/records/18907380
Please use the above DOI when citing the specific version of the code used in your research.
