#!/usr/bin/env bash

#framerate" how long each image is shown in seconds 
#start_number: the first image to index
#-i: jpg images inputed 
#-c:v video encoder used
#crf : quality vs file size
#pix_fmt: pixel format

ffmpeg -framerate 1\
  -start_number 0\
  -i "/home/casakurai/scratch/asyn-phase-binding-data/Figures/MD-images/flat/defect-cut-off-10A/frame%d.png" \
  -vf "scale=ceil(iw/2)*2:ceil(ih/2)*2,format=rgb24" \
  -c:v libx264  \
  -crf 20 \
  -pix_fmt yuv420p \
  "/home/casakurai/scratch/asyn-phase-binding-data/Figures/MD-images/flat/compiled/flat-compiled-defect-data-cutoff-10A.mp4"

# Delete the frames after ffmpeg finishes
#rm final.top-DPPC-DPPA-strain.2-lipidtype.*.jpg
