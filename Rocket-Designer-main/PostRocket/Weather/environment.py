# defines the environment class, which contains all of the weather data for a given time and place

import time as tm
from matplotlib import pyplot as plt
import numpy as np
import fetch

# Constants
heightSteps: list[int] = [
    110,
    320,
    500,
    800,
    1000,
    1500,
    1900,
    3000,
    4200,
    5600,
    7200,
    9600,
]

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
	
	def getIndices(self, height) -> tuple[int, int]:
		for i in range(len(heightSteps)):
			if heightSteps[i] > height:
				return i-1, i
	
	def getAtHeight(self, property: str, height: int, hour: int):
		fetch.validateProperties([property])
		if not 0 <= height <= heightSteps[-1]:
			raise Exception(f"Height of {height} out of bounds")
		if not 0 <= hour <= self.atmosphere[property].shape[0]:
			raise Exception(f"Hour {hour} is out of bounds\nMax hour for this API call is {self.atmosphere[property].shape[0]}")
		
		if height < heightSteps[0]:
			return self.atmosphere[property][hour,0]
		
		if height in heightSteps:
			return self.atmosphere[property][hour,heightSteps.index(height)]
		
		lower, higher = self.getIndices(height)

		# use those indices to calculate linear approximation between the heights at those indices

		raise NotImplementedError


# Just for testing, probably remove
def main():
	# env = Environment(32.938358, -106.912406)
	# env.fetch_openMeteoData(["temp"], 2)
	# print(env.atmosphere["windSpeed"][23,:])
	plt.plot(list(range(1,13)), heightSteps)
	# plt.plot(np.arange(1,13), 60 * np.arange(1,13)**2.2)
	plt.show()
	# print(np.array([1,2,3])**2)

if __name__ == "__main__":
	main()
