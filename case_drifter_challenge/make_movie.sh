#! /bin/bash

ffmpeg -framerate 25 -pattern_type glob  -i "images_hindcast/map_*.png" output.mp4
