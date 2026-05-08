#!/usr/bin/env bash

#framerate" how long each image is shown in seconds 
#start_number: the first image to index
#-i: jpg images inputed 
#-c:v video encoder used
#crf : quality vs file size
#pix_fmt: pixel format
#  -vf "scale=ceil(iw/2)*2:ceil(ih/2)*2,format=rgb24" \
# for the MD simulations the format needs to be .tga
movietype="bead-type"
system="DOPC"
configuration="flat"
view="side"

ffmpeg -framerate 1\
  -start_number 0\
  -i "/home/casakurai/scratch/asyn-phase-binding-data/Figures/MD-images/$movietype/$configuration/$system/$view/frame%04d.tga" \
  -c:v libx264  \
  -crf 20 \
  -pix_fmt yuv420p \
  "/home/casakurai/scratch/asyn-phase-binding-data/Figures/MD-images/lipid-type/$configuration/compiled/individual-systems/$system-$configuration-$view-MD-lipid-distribution.mp4"

# Delete the frames after ffmpeg finishes
#rm final.top-DPPC-DPPA-strain.2-lipidtype.*.jpg

