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

		self.temperature: np.ndarray = np.zeros(1);
		self.pressure: np.ndarray = np.zeros(1);
		self.density: np.ndarray = np.zeros(1);
		self.windSpeed: np.ndarray = np.zeros(1);
		self.windDirection: np.ndarray = np.zeros(1);
		self.gravity: np.ndarray = np.zeros(1);

	def calculate_gravity(self):
		raise NotImplementedError
	
	def fetch_openMeteoData(self):
		temp = fetch.fetch_data(self.lat, self.lon, "current")


# Just for testing, probably remove
def main():
	print()

if __name__ == "__main__":
	main()
