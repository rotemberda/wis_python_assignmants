from math import pi

#asks for radius variable and validate it
while True:
    try:
        radius = int(input("What's the radius? "))
        break
    except:
        print("Invalid input! Input type isn't 'int'")

#print the aera and circumference of the circle

print("The area of the circle is", radius ** 2 * pi)
print("The circumference of the circle is", radius * 2 * pi)