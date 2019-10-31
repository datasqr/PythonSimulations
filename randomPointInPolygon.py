import numpy
from random import choices
import random
import csv
import math
import pandas as pd 
from datetime import datetime, timedelta
import names
import json
import geog
import shapely.geometry
from shapely.geometry import Polygon, Point, LineString
import multiprocessing
import time

# Function to create polygon (circle)
# Input: coordinates of the polygon center (lat, lng) and its radious
# Output: Polygon in a form of circle
def createPolygon(lat, lng, radious):
    p = shapely.geometry.Point([lat, lng])
    n_points = 20
    d = radious  # meters
    angles = numpy.linspace(0, 360, n_points)
    polygon = geog.propagate(p, angles, d)
    result = shapely.geometry.mapping(shapely.geometry.Polygon(polygon))
    return Polygon(result["coordinates"][0])

# Function to get random location of the point inside the polygon
# Input: Polygon
# Outpuy: Cooridination of the random point
def random_point_within_Polygon(poly):
    min_x, min_y, max_x, max_y = poly.bounds

    x = random.uniform(min_x, max_x)
    x_line = LineString([(x, min_y), (x, max_y)])
    x_line_intercept_min, x_line_intercept_max = x_line.intersection(poly).xy[1].tolist()
    y = random.uniform(x_line_intercept_min, x_line_intercept_max)

    return Point([x, y])

polygon = createPolygon(21.01, 52.22, 5000)
feature1 = [{'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(polygon)}]
point = random_point_within_Polygon(polygon)
feature2 = [{'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(point)}]

with open('somefile1.txt', 'a') as the_file:
    the_file.write(str(feature1))

with open('somefile2.txt', 'a') as the_file:
    the_file.write(str(feature2))
