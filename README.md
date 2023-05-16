# CSreconstruction

This repository provides code for the compressed sensing reconstruction of single channel MRI data. It was developed in the context of preclinical fluorine-19 MRI and provides both very fast computation and automatic regularization strength optimization. 

The implemented algorithm is Goldstein et al.'s accelerated alternating direction method of multipliers. Very fast reconstruction with quick and reliable covergence is achieved by using a Fourier transform-based exact matrix inversion where other algorithms rely on conjugate gradient optimization. This exploits the circulant matrix structure arising from single channel data with periodic boundary conditions.

In the automatic version, multiple reconstructions are performed optimizing the regularization strength until a desired deviation of the reconstruction from the measured data is achieved (discrepancy principle).

For details, please see 

    Starke, Ludger, et al. 
    "Performance of compressed sensing for fluorine‐19 magnetic resonance 
    imaging at low signal‐to‐noise ratio conditions." 
    Magnetic resonance in medicine 84.2 (2020): 592-608.

I also included the relevant parts of my master thesis deriving the used update equations. The most relevant sections are 3.1.2, 3.1.3 and A.2.

This code was developed back in 2017/18, but as the algorithm is still effective and competetitive for its intended application, I wanted to make it accessible.

The 3 scripts (exampleCSreco_2D.m, exampleCSreco_19F.m, exampleCSreco_1H.m) are identical except for the extracted data and some default parameter values.

If you use this software, please cite the above publication. The code is licensed under GNU GPLv3.
