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

coordinates: dict[str, tuple[float, float]] = {
    "Spaceport America": (32.938358, -106.912406),
    "Utah1": (37.931728, -113.053677),
    "Utah2": (37.945524, -113.033278),
    "Texas": (31.049802, -103.547313),
}


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

    def fetchOpenMeteoData(
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
            raise Exception(
                f"{property} is valid, but was not given as a property to be retrieved from the API"
            )
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


# Mostly for testing, but it does show the proper use of some of the functions so I guess I'll leave it
def main():
    env = Environment(*coordinates["Texas"])
    env.fetchOpenMeteoData(fetch.validProperties, days=1)
    print(env.getAtHeight("temp", 750, 1))


if __name__ == "__main__":
    main()
