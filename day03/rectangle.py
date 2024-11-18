#asks input from user and check if valid
while True:
    try:
        width = int(input("what the width? "))
        break
    except:
        print("Invalid input! Input type isn't 'int'")

while True:
    try:
        length = int(input("what the length? "))
        break
    except:
        print("Invalid input! Input type isn't 'int'")

#print the aera and circumference of the rectangle

print("The area of the rectangle is", length * width)
print("The circumference of the rectangle is", length * 2 + width * 2)