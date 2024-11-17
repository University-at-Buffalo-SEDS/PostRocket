# defines the environment class, which contains all of the weather data for a given time and place

import time as tm
import numpy as np
import fetch

# Constants
# idk what constants we need but they'll go here

class Environment:
	
	"""
	Initialize weather data given a specific latitude, longitude, and time.
	Currently tracks temperature, pressure, density, wind direction, and
	wind magnitude at stepped altitudes

	!!These variables will probably be ndarrays but these are obviously just placeholders at the moment
	"""
	def __init__(self, lat: float, lon: float) -> None:
		self.lat = lat
		self.lon = lon
		self.atmosphere: dict[str, np.ndarray]

	def calculate_gravity(self):
		raise NotImplementedError
	
	# properties should be from list ["temp", "humidity", "windSpeed", "windDirection"].
	# Days is number of days in the futureforecast
	def fetch_openMeteoData(self, properties: list[str], days: int) -> dict[str, np.ndarray]:
		self.atmosphere = fetch.fetch_data(self.lat, self.lon, properties, days)
		return self.atmosphere


# Just for testing, probably remove
def main():
	env = Environment(32.938358, -106.912406)
	env.fetch_openMeteoData(["temp", "humidity", "windSpeed", "windDirection"], 1)
	print(env.atmosphere["windSpeed"][23,:])

if __name__ == "__main__":
	main()
