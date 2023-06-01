# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:48:20 2023

@author: simpson-c
"""

############################# Running Starch4Kranz:
# Go to directory
import os
import matplotlib.pyplot as plt

# Running Starch4Kranz
# Run this in the same directory as images
os.chdir(r"D:/Conor-temp/Maize")

image_name = []
vd = []
processing = []

for filename in os.listdir('.'):
    if filename.endswith('.tif'):
        vein_density, Processed, image =  Starch4Kranz(filename = filename, split_mode = "monocot", trim_factor = 90, pixel_length_um =0.7, x_pixels = 1536, y_pixels = 2048, initial_blur_factor=50)
        
        image_name.append(filename)
        vd.append(vein_density)
        processing.append(Processed)
        
        plt.imshow(image)
        plt.title(Processed)
        plt.pause(0.1)  # Pause to display the figure briefly before moving to the next iteration

plt.show()


results = [['Image_Name', 'VD_um_mm2', 'Processing']]
results += list(zip(image_name, vd, processing))

header = results[0]
results = results[1:]

results = [header] + sorted(results)



