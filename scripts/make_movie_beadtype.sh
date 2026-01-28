#!/usr/bin/env bash

#framerate" how long each image is shown in seconds 
#start_number: the first image to index
#-i: jpg images inputed 
#-c:v video encoder used
#crf : quality vs file size
#pix_fmt: pixel format
#  -vf "scale=ceil(iw/2)*2:ceil(ih/2)*2,format=rgb24" \
# for the MD simulations the format needs to be .tga
movietype="lipid-type"
configuration="flat"
view="side"

ffmpeg -framerate 1\
  -start_number 0\
  -i "/home/casakurai/scratch/asyn-phase-binding-data/Figures/MD-images/$movietype/$configuration/compiled/all-systems/$view-frame%d.png" \
  -c:v libx264  \
  -vf "scale=ceil(iw/2)*2:ceil(ih/2)*2,format=rgb24" \
  -crf 20 \
  -pix_fmt yuv420p \
  "/home/casakurai/scratch/asyn-phase-binding-data/Figures/MD-images/$movietype/$configuration/compiled/all-systems/$view-MD-$movietype.mp4"

# Delete the frames after ffmpeg finishes
#rm final.top-DPPC-DPPA-strain.2-lipidtype.*.jpg

