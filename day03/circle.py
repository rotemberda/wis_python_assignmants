from math import pi
import argparse

#asks for radius variable and validate it
parser = argparse.ArgumentParser()
parser.add_argument('--radius', help= 'Radius of circle', required= True, type = int)

args = parser.parse_args()

radius = args.radius

#print the aera and circumference of the circle

print("The area of the circle is", radius ** 2 * pi)
print("The circumference of the circle is", radius * 2 * pi)