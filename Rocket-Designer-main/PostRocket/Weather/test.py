#The following code was NOT MADE BY ME (so don't expect me to explain how it works), a modified version is found at "https://open-meteo.com/en/docs#current=temperature_850hPa,temperature_800hPa,temperature_700hPa,wind_speed_850hPa,wind_speed_800hPa,wind_speed_700hPa,wind_direction_850hPa,wind_direction_800hPa,wind_direction_700hPa&temperature_unit=fahrenheit&wind_speed_unit=ms&precipitation_unit=inch&timezone=America%2FDenver&forecast_days=1"

from datetime import datetime
import openmeteo_requests

import requests_cache
import pandas as pd
from retry_requests import retry

# Setup the Open-Meteo API client with cache and retry on error
cache_session = requests_cache.CachedSession('.cache', expire_after = 3600)
retry_session = retry(cache_session, retries = 5, backoff_factor = 0.2)
openmeteo = openmeteo_requests.Client(session = retry_session)

latitude = 32.938358
longitude = -106.912406

# Make sure all required weather variables are listed here
# The order of variables in hourly or daily is important to assign them correctly below
url = "https://api.open-meteo.com/v1/forecast"
params = {
	"latitude": 32.938358,
	"longitude": -106.912406,
	"current": ["temperature_850hPa", "temperature_800hPa", "temperature_700hPa", "wind_speed_850hPa", "wind_speed_800hPa", "wind_speed_700hPa", "wind_direction_850hPa", "wind_direction_800hPa", "wind_direction_700hPa"],
	"temperature_unit": "fahrenheit",
	"wind_speed_unit": "ms",
	"precipitation_unit": "inch",
	"timezone": "America/Denver",
	"forecast_days": 1
}
responses = openmeteo.weather_api(url, params=params)

# Process first location. Add a for-loop for multiple locations or weather models
response = responses[0]
print(f"Coordinates {response.Latitude()}°N {response.Longitude()}°E")
print(f"Elevation {response.Elevation()} m asl")
print(f"Timezone {response.Timezone()} {response.TimezoneAbbreviation()}")
print(f"Timezone difference to GMT+0 {response.UtcOffsetSeconds()} s")

# Process current data. The order of variables needs to be the same as requested.
current = response.Current()
current_temperature_850hPa = current.Variables(0).Value()
current_temperature_800hPa = current.Variables(1).Value()
current_temperature_700hPa = current.Variables(2).Value()
current_wind_speed_850hPa = current.Variables(3).Value()
current_wind_speed_800hPa = current.Variables(4).Value()
current_wind_speed_700hPa = current.Variables(5).Value()
current_wind_direction_850hPa = current.Variables(6).Value()
current_wind_direction_800hPa = current.Variables(7).Value()
current_wind_direction_700hPa = current.Variables(8).Value()

#"https://api.open-meteo.com/v1/forecast?latitude=52.52&longitude=13.41&current=temperature_2m,relative_humidity_2m,apparent_temperature,is_day&temperature_unit=fahrenheit&wind_speed_unit=ms&precipitation_unit=inch&timezone=America%2FDenver&forecast_days=1"

time_date = datetime.fromtimestamp(current.Time())

print(f"Current time {time_date}")
print(f"Current temperature_850hPa {current_temperature_850hPa}")
print(f"Current temperature_800hPa {current_temperature_800hPa}")
print(f"Current temperature_700hPa {current_temperature_700hPa}")
print(f"Current wind_speed_850hPa {current_wind_speed_850hPa}")
print(f"Current wind_speed_800hPa {current_wind_speed_800hPa}")
print(f"Current wind_speed_700hPa {current_wind_speed_700hPa}")
print(f"Current wind_direction_850hPa {current_wind_direction_850hPa}")
print(f"Current wind_direction_800hPa {current_wind_direction_800hPa}")
print(f"Current wind_direction_700hPa {current_wind_direction_700hPa}")
