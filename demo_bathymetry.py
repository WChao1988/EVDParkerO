from faatobathymetry import EVDBATHYMETRY
import numpy as np





if __name__ == "__main__":
    lat_up, lat_down, lon_left, lon_right = 19, 16, 130, 133
    # lat_up, lat_down, lon_left, lon_right = 18, 15, -164, -161
    # lat_up, lat_down, lon_left, lon_right = 36, 33, 143, 146
    gravity = np.loadtxt(f'./data_require_bathymetry_inversion/gravity_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    water = np.loadtxt(f'./data_require_bathymetry_inversion/water_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    sediment = np.loadtxt(f'./data_require_bathymetry_inversion/sediment_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    longwave_topo = np.loadtxt(f'./data_require_bathymetry_inversion/longwave_topo_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    longrkm = 333.5847799336763 #333.5847799336763 #333.584779933676
    longckm = 318.1421721146795 #319.8447222955345 # 274.9058765744849

    parameters = {
        "faa_matrix": gravity,
        "sediment_matrix": sediment,
        "delta_rhos": {0: 1.3, 1: 0.85},
        "delta_rho_initial": 1.77,
        "mus": {0: 0.01, 1: 0.01},
        "reference_depths": {0: -(water-sediment).mean(), 1: -water.mean()},
        "longrkm": longrkm,
        "longckm": longckm,
        "wh": 0.2,
        "alpha": 8,
        "longwave_topo": longwave_topo
    }
    BATHYMETRY = EVDBATHYMETRY(**parameters)
    bathymetry = BATHYMETRY.downward(t=6)

