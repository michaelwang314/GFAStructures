# GFAStructures

<p align="center">
    <img src="https://github.com/michaelwang314/GFAStructures/blob/main/examples/structures_example.png" width="750">
</p>

## Introduction

The purpose of this package is to simulate the relaxation of geometrically frustrated assemblies made up of misfitting subunits.  One can define the shape of the subunits (via the positions of their interaction sites) and their inter-subunit interactions (e.g. Harmonic springs, Lennard-Jones) with certain connectivities. See ... for more details as well as example codes in the examples folder.

Setting up the simulation is quite simple:
1) Create a subunit as a rigid body composed of a number of interaction sites
2) Define the pairs of interacting sites via an interaction matrix
3) Generate an initial arrangement of subunits with given connectivities (e.g. in a sheet or on a cylinder)
4) Define the interactions (e.g. Harmonic springs, Lennard-Jones) between pairs of interacting sites
5) Run simulation and relax the structure

## Downloading package

To download, open Julia command line, press `]`, and use `add`
```
(v1.10) pkg> add https://github.com/michaelwang314/GFAStructures
```
or do
```
julia> using Pkg
julia> Pkg.add(url = "https://github.com/michaelwang314/GFAStructures")
```

## Creating subunits

...TBC