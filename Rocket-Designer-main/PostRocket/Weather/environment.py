# defines the environment class, which contains all of the weather data for a given time and place

from logging import warning
from sys import exception
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
    Currently tracks temperature, pressure, wind direction, and
    wind speed at stepped altitudes
    """

    def __init__(self, lat: float, lon: float) -> None:
        self.lat = lat
        self.lon = lon
        self.atmosphere: dict[str, np.ndarray] = {}

    def fetch_openMeteoData(
        self, properties: list[str], days: int
    ) -> dict[str, np.ndarray]:
        """
        Calls the openMeteo API for the desired properties
        Properties should be from list ["temp", "humidity", "windSpeed", "windDirection"]
        Days is number of days in the future to forecast, and should be in the list
        """
        self.atmosphere = fetch.fetch_data(self.lat, self.lon, properties, days)
        return self.atmosphere

    def getIndices(self, height) -> tuple[int, int]:
        """Get the indices of the height steps below and above the given height"""
        for i in range(len(heightSteps)):
            if heightSteps[i] > height:
                return i - 1, i

    def getAtHeight(self, property: str, height: int, hour: int):
        """
        Return a given property from the atmosphere attribute given a specific height and time
        If height is below the minimum height step, returns the value at the minimum height step
        Otherwise, linearly approximates the property between the two nearest height steps
        """
        if self.atmosphere == {}:
            raise Exception("Call the API before accessing environment data")
        fetch.validateProperties([property])
        if property not in self.atmosphere.keys():
            raise Exception(f"{property} is valid, but was not given as a property to be retrieved from the API")
        if not 0 <= height <= heightSteps[-1]:
            raise Exception(f"Height of {height} out of bounds")
        if not 0 <= hour < self.atmosphere[property].shape[0]:
            raise Exception(
                f"Hour {hour} is out of bounds\nMax hour for this API call is {self.atmosphere[property].shape[0]}"
            )

        if height < heightSteps[0]:
            return self.atmosphere[property][hour, 0]

        if height in heightSteps:
            return self.atmosphere[property][hour, heightSteps.index(height)]

        lowIdx, highIdx = self.getIndices(height)
        low = heightSteps[lowIdx]
        high = heightSteps[highIdx]

        # use those indices to calculate linear approximation between the heights at those indices
        return (
            self.atmosphere[property][hour, highIdx]
            - self.atmosphere[property][hour, lowIdx]
        ) * ((height - low) / (high - low)) + self.atmosphere[property][hour, lowIdx]


# Just for testing, probably remove
def main():
    # env = Environment(32.938358, -106.912406)
    # env.fetch_openMeteoData(["t", "temp", "pawfw", "pressure", "humidity"], 1)
    # print(env.atmosphere["windSpeed"][23,:])
    # plt.plot(list(range(1, 13)), heightSteps)
    # plt.plot(np.arange(1,13), 60 * np.arange(1,13)**2.2)
    # plt.show()
    # print(np.array([1,2,3])**2)
    # env = Environment(32.938358, -106.912406)
    # env.atmosphere["temp"] = np.array([[x*100 for x in range(1,13)]])
    # print(env.getAtHeight("temp", 9600, 1))


if __name__ == "__main__":
    main()
