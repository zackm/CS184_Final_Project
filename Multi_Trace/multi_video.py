#!/usr/bin/env python

import os

# go through each obj file in input folder and raytrace it

# array of scene numbers
scenes = [1,2,3,4,5]

for num in scenes:
    folder_name = str(num) + "_images"
    os.system("../final -d -output "+folder_name+"_point")
    os.system("../final -i -output "+folder_name+"_part")
    os.system("../final -p -output "+folder_name+"_iso")

