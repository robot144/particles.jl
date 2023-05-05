
# Introduction

This folder contains scripts for simulation of several drifters that were deployed in March and October 2017 in the German Bight. The dataset originates from:

- J. Meyerjurgens, T. H. Badewien, S. P. Garaba, J.-O. Wol, and O. Zielinski, A state-of-the-art compact surface drifter reveals pathways of oating marine litter in the german bight," Frontiers in
Marine Science, vol. 6, p. 58, 2019.

## How to run

The standard runs in the script require access to flow fields from the 3D DCSM-FM model. These are too large (10x 30Gb) for github, so they are stored elsewhere. They can be downloaded with the 
script 'download_data.sh' using bash and wget on linux or 'download_data.jl' using julia.

## Ongoing work

Add 3d interpolation for delft3d flexible mesh.
- [ ] add unit tests for conversion of delft3dfm to zarr
- [ ] add lock exchange map file to unit test
- [ ] add interpolation with FullGrid

