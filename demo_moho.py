from faatomoho import EVDMOHO
import numpy as np




if __name__ == "__main__":
    lat_up, lat_down, lon_left, lon_right = 20, 15, 129, 134
    # lat_up, lat_down, lon_left, lon_right = 19, 14, -165, -160
    # lat_up, lat_down, lon_left, lon_right = 38, 31, 141, 148
    gravity = np.loadtxt(f'./data_require_moho_inversion/gravity_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    water = np.loadtxt(f'./data_require_moho_inversion/water_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    sediment = np.loadtxt(f'./data_require_moho_inversion/sediment_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    crust = np.loadtxt(f'./data_require_moho_inversion/crust_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    moho = np.loadtxt(f'./data_require_moho_inversion/moho_{lat_up}_{lat_down}_{lon_left}_{lon_right}.txt')
    longrkm = 555.9746332227937 #555.9746332227937  778.3644865119115
    longckm = 530.2272075400947 #533.0657964209347  641.3424465308995


    parameters = {
        'faa_matrix': gravity - gravity.mean(),
        'delta_bnds': {0: None, 1: (water - sediment) - (water - sediment).mean()},
        'delta_rhos': {0: 2.5, 1: 1.63},
        'mus': {0: 0, 1: 0.01},
        'reference_depths': {0: -moho.mean(), 1: -crust.mean()},
        'longrkm': longrkm,
        'longckm': longckm,
        'wh': 0.2,
        "alpha": 64
    }

    MOHO = EVDMOHO(**parameters)
    moho_predicted = MOHO.downward(t=4)


