"""
Fetches data from APIs
Currently (and probably permanently) only uses OpenMeteo API
"""

import environment
import openmeteo_requests
import requests_cache
from retry_requests import retry
import numpy as np
from openmeteo_sdk.WeatherApiResponse import WeatherApiResponse, VariablesWithTime

# Used to determine a valid or invalid call to the fetchData method. Also used from environment.py
validProperties: list[str] = ["temp", "humidity", "windSpeed", "windDirection"]
validDays: list[int] = [1, 3, 7, 14, 16]

# Used in construction of the API call
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


def validateProperties(properties: list[str]) -> None:
    """Determines is given properties are valid, throw exception if not. This isn't very good practice"""
    validSet = set(validProperties)
    actualSet = set(properties)
    if not actualSet.issubset(validSet):
        raise Exception(
            f"{list(actualSet - validSet)} are not valid properties\nValid properties are {validProperties}"
        )


def validateDays(days: int) -> None:
    """Determines if given number if days is valid, throws exception if not. This isn't very good practice"""
    if days not in validDays:
        raise Exception(
            f"{days} is not a valid forecast range\nDays must be 1, 3, 7, 14, or 16"
        )


# Calls the open-meteo api and returns a dictionary of each of the requested properties
def fetch_data(
    lat: float, lon: float, properties: list[str], days: int
) -> dict[str, np.ndarray]:
    """Fetch data from openMeteo API"""

    # This isn't very good practice and should probably be changed
    validateDays(days)
    validateProperties(properties)

    cache_session = requests_cache.CachedSession(".cache", expire_after=3600)
    retry_session = retry(cache_session, retries=5, backoff_factor=0.2)
    openmeteo = openmeteo_requests.Client(session=retry_session)

    # Construct the properties that are desired in the necessary format for the API call
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

    hourly: VariablesWithTime = response.Hourly()

    numHeightSteps = len(environment.heightSteps)
    hours = days * 24

    out: dict[str, np.ndarray] = {}
    for property in properties:
        out[property] = np.zeros(
            (hours, numHeightSteps)
        )  # Each property in the output will be (hours x heightSteps) size

    for i in range(len(properties)):  # for each property
        var = np.zeros((hours, numHeightSteps))
        for j in range(
            numHeightSteps * i, numHeightSteps * (i + 1)
        ):  # for each height step of that property
            var[:, j % numHeightSteps] = hourly.Variables(j).ValuesAsNumpy()
        out[properties[i]] = var

    return out
