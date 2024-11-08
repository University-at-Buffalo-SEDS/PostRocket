"""
Fetches data from APIs
Currently (and probably permanently) only uses OpenMeteo API
"""

from datetime import datetime
import openmeteo_requests
import requests_cache
import pandas as pd
from retry_requests import retry
import numpy as np
import requests
import time as tm


# I'm not totally sure about the return type, or whether data analysis will be
# done in this file or the environment class. Probably the latter
def fetch_data(lat: float, lon: float, timeType: str, ) -> np.ndarray:
    cache_session = requests_cache.CachedSession(".cache", expire_after=3600)
    retry_session = retry(cache_session, retries=5, backoff_factor=0.2)
    openmeteo = openmeteo_requests.Client(session=retry_session)

    url = "https://api.open-meteo.com/v1/forecast"
    


    params = {
        "latitude": lat,
        "longitude": lon,
        timeType: [
            "temperature_850hPa",
            "temperature_800hPa",
            "temperature_700hPa",
            "wind_speed_850hPa",
            "wind_speed_800hPa",
            "wind_speed_700hPa",
            "wind_direction_850hPa",
            "wind_direction_800hPa",
            "wind_direction_700hPa",
        ],
        "temperature_unit": "fahrenheit",
        "wind_speed_unit": "ms",
        "timezone": "America/Denver",
        "forecast_days": 1,
    }

    return np.zeros(1)
