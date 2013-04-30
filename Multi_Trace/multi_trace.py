#!/usr/bin/env python

import os

# go through each obj file in input folder and raytrace it

# array of filenames
filenames = [name for name in os.listdir('./input_obj/')]
for i in range(0,len(filenames)):
    print "Working on: " + filenames[i]
    name = filenames[i]
    name_pic = name[0:-3] + "png"
    os.system("./start_trace "+name+" "name_pic+" 1000 1000")
