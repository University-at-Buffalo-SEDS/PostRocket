"""
Fetches data from APIs
"""

import numpy as np
import requests
import time as tm

# I'm not totally sure about the return type, or whether data analysis will be 
# done in this file or the environment class. Probably the latter
def fetch_data(lat: float, lon: float, time: tm.struct_time) -> np.ndarray:
	return np.zeros(1)


