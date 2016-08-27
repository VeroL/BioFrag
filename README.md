# BioFrag | Edge response
The BioFrag Edge response software is a collection of novel image processing and computational methods designed for Landscape Ecologists, to measure habitat framgentation and analyse its impact on species abundance and biodiversity from species census data and remotely sensed images.

In an endeavour to improve on current methods used to assess habitat fragmentation we have developed a continuous
 spatial measure of Edge Influence (EI), which we implemented in this software. 
Contrary to the distance to nearest edge measure, the design of our Edge Influence measure takes into account important 
but often overlooked characteristics of fragmented landscapes: distance to multiple edges and edge shape, fragment contrast
 and habitat quality. Additionally, EI also accounts for a set of commonly used metrics: fragment size, fragment isolation 
and local configuration. Our EI measure may be used as a local proxy for fragmentation instead of a plurality of interdependent metrics
 as it encompasses edge effects, local habitat amount, local contrast, shape and local connectivity. 
The Edge_response software computes a map of EI from a map of a continous variable representing land cover (e.g. % Tree Cover).


Other novel tools include:

-automated classification of species abundance patterns in response to EI and land cover

-mapping of predicted species abundances from EI and land cover

-computation of species response metrics: impact of fragmentation on species abundance, and species sensitivity to EI

-computation of landscape metrics: configuration and amount of edge effects

Edge_response software is implemented in Matlab 2012a and distributed with the Matlab Runtime environment.



24/08/2016 

Spatial methods developement and implementation by Dr Véronique Lefebvre, Dr Marion Pfeifer and Dr. Robert Ewers

Forest Ecology and Conservation group : http://forestecology.net/

Imperial College London
