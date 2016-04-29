Spline based smoothing, 2016-01-25
===========================

This directory contains all of the analysis files required to reproduce the figures and results in the manuscript.

### bspline.m

This is the basic route for creating b-splines of arbitrary order, at arbitrary knot points, on an arbitrary grid. It is based off of the algorithm in *A Practical Guide to Splines* by Carl de Boor.

### LoadFigureDefaults.m

A file with the default figure size settings that is referenced by all MakeFigure scripts.

# Manuscript Figures and Tables

### MakeBSplineFigure.m

Creates a figure to show an example b-spline and its derivatives at a few orders.

### MakeInterpolationFigure.m

Creates a figure to show how b-splines can be used to interpolate data.

### MakeOptimalParameterFigure.m

Strides through synthetic data created by MakeAndSaveMaternNoise.m to find the optimal tension parameter, and then compares to our two methods.

