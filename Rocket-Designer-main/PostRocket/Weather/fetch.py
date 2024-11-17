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
from openmeteo_sdk.WeatherApiResponse import WeatherApiResponse, VariablesWithTime

propertiesMap: dict[str, list[str]] = {
    "temp": [
        "temperature_1000hPa",  # 110 m
        "temperature_975hPa",  # 320 m
        "temperature_950hPa",  # 500 m
        "temperature_925hPa",  # 800 m
        "temperature_900hPa",  # 1000 m
        "temperature_850hPa",  # 1500 m
        "temperature_800hPa",  # 1900 m
        "temperature_700hPa",  # 3000 m
        "temperature_600hPa",  # 4200 m
        "temperature_500hPa",  # 5600 m
        "temperature_400hPa",  # 7200 m
        "temperature_300hPa",  # 9600 m
    ],
    "humidity": [
        "relative_humidity_1000hPa",  # 110 m
        "relative_humidity_975hPa",  # 320 m
        "relative_humidity_950hPa",  # 500 m
        "relative_humidity_925hPa",  # 800 m
        "relative_humidity_900hPa",  # 1000 m
        "relative_humidity_850hPa",  # 1500 m
        "relative_humidity_800hPa",  # 1900 m
        "relative_humidity_700hPa",  # 3000 m
        "relative_humidity_600hPa",  # 4200 m
        "relative_humidity_500hPa",  # 5600 m
        "relative_humidity_400hPa",  # 7200 m
        "relative_humidity_300hPa",  # 9600 m
    ],
    "windSpeed": [
        "wind_speed_1000hPa",  # 110 m
        "wind_speed_975hPa",  # 320 m
        "wind_speed_950hPa",  # 500 m
        "wind_speed_925hPa",  # 800 m
        "wind_speed_900hPa",  # 1000 m
        "wind_speed_850hPa",  # 1500 m
        "wind_speed_800hPa",  # 1900 m
        "wind_speed_700hPa",  # 3000 m
        "wind_speed_600hPa",  # 4200 m
        "wind_speed_500hPa",  # 5600 m
        "wind_speed_400hPa",  # 7200 m
        "wind_speed_300hPa",  # 9600 m
    ],
    "windDirection": [
        "wind_direction_1000hPa",  # 110 m
        "wind_direction_975hPa",  # 320 m
        "wind_direction_950hPa",  # 500 m
        "wind_direction_925hPa",  # 800 m
        "wind_direction_900hPa",  # 1000 m
        "wind_direction_850hPa",  # 1500 m
        "wind_direction_800hPa",  # 1900 m
        "wind_direction_700hPa",  # 3000 m
        "wind_direction_600hPa",  # 4200 m
        "wind_direction_500hPa",  # 5600 m
        "wind_direction_400hPa",  # 7200 m
        "wind_direction_300hPa",  # 9600 m
    ],
}


# Calls the open-meteo api and returns a dictionary of each of the requested properties
def fetch_data(
    lat: float, lon: float, properties: list[str], days: int
) -> dict[str, np.ndarray]:
    
    cache_session = requests_cache.CachedSession(".cache", expire_after=3600)
    retry_session = retry(cache_session, retries=5, backoff_factor=0.2)
    openmeteo = openmeteo_requests.Client(session=retry_session)

    callProperties = [propertiesMap[prop] for prop in properties]

    url = "https://api.open-meteo.com/v1/forecast"

    params = {
        "latitude": lat,
        "longitude": lon,
        "hourly": callProperties,
        "temperature_unit": "fahrenheit",
        "wind_speed_unit": "ms",
        "timezone": "America/Denver",  # I don't know if this is necessary but it currently doesn't change if we change the coords
        "forecast_days": days,
    }

    response: WeatherApiResponse = openmeteo.weather_api(url, params=params)[0]
    # print(f"Coordinates {response.Latitude()}°N {response.Longitude()}°E")
    # print(f"Elevation {response.Elevation()} m asl")

    hourly: VariablesWithTime = response.Hourly()
    # print(hourly.Variables(0).ValuesAsNumpy())

    heightSteps = 12
    hours = days * 24

    out: dict[str, np.ndarray] = {}
    for property in properties:
        out[property] = np.zeros((hours, heightSteps)) # Each property in the output will be (hours x heightSteps) size


    for i in range(len(properties)): # for each property
        var = np.zeros((hours, heightSteps)) # for each height step of that property
        for j in range(heightSteps * i, heightSteps * (i + 1)):
            var[:, j % heightSteps] = hourly.Variables(j).ValuesAsNumpy()
        out[properties[i]] = var

    return out
