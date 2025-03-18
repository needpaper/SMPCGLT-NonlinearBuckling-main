import os
import csv
import math

def createside(Sectionalradius, angle, radius, Lumbus):
    sr = Sectionalradius
    r1 = radius
    a1 = angle
    w = Lumbus
    r2 = (sr - r1 + r1 * math.cos(angle * math.pi / 180.0)) / (1 - math.cos(angle * math.pi / 180.0))
    Design_space = 120.0
    w1 = Design_space * 0.5 - angle * r1 * math.pi / 180 - angle * r2 * math.pi / 90 - w * 0.5
    return r2, w1

def find_valid_combinations():
    valid_combinations = []
    for Sectionalradius in range(12, 19):
        for angle in range(30, 81, 5):  
            for radius in range(10, Sectionalradius + 1):
                for Lumbus in range(1, 11):
                    r2, w1 = createside(Sectionalradius, angle, radius, Lumbus)
                    if r2 > 2 and w1 > 2:
                        valid_combinations.append((Sectionalradius, angle, radius, Lumbus))
    return valid_combinations

def save_combinations_to_csv(combinations, filename='valid_combinations.csv'):
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(['Sectionalradius', 'angle', 'radius', 'Lumbus'])
        # Write data
        writer.writerows(combinations)

if __name__ == "__main__":
    valid_combinations = find_valid_combinations()
    save_combinations_to_csv(valid_combinations)