[![DOI](https://zenodo.org/badge/279592337.svg)](https://zenodo.org/badge/latestdoi/279592337)
# The calculation of the Hilbert transform by a non-equispaced Fast Fourier Transform (NFFT) and its application to gravity and magnetic data
Computes the Hilbert transform of vector fields in the space domain (see Nabighian, 1984). The convolution sum is efficiently calculated by (Potts et. al., 2004)s' fastsum algorithm.
 This package includes the following subdirectories:
 - 1DHilbertConvolution: Computes the 1D Hilbert transform with two different convolution approaches. Uses a gravimetric field as an example. 
 - 1DHilbertFastsum: Computes the 1D Hilbert transfrom with an efficient fastsum algorithm. Uses a gravimetric field as an example.
 - 3DHilbertFastsum: Computes 3D Hilbert transfrom with an efficient fastsum algorithm. Uses a magnetic field as an example.
 - 3DHilbertFastsumReal: Computes 3D Hilbert transfrom with an efficient fastsum algorithm. Uses uses real world magnetic data as an example.
 
 In addition, the frequency approach for calculating the Hilbert transform is used for comparison in all the subdirectories mentioned above.
 
# Third party packages contained in this repository:
- Tearle, Matt (2020). Software: Simple real Fourier series approximation. https://www.mathworks.com/matlabcentral/fileexchange/31013-simple-real-fourierseries-approximation. Downloaded: 20 February 2020.
- Semechko, Anton (2019). Software: Bounding Spheres and Circles. https://github.com/AntonSemechko/Bounding-Spheres-And-Circles. Downloaded: 5 November 2019.
- Potts, Daniel and Stefan Kunis (2019). Software: Nonequispaced fast Fourier transform - Matlab/Octave interface. https://www-user.tu-chemnitz.de/~potts/nfft/download.php. Downloaded: 9 August 2019.
- Stafford, Roger (2005). Software: Random Points in an n-Dimensional Hypersphere. https://www.mathworks.com/matlabcentral/fileexchange/9443-random-points-in-an-n-dimensional-hypersphere. Downloaded: 30 October 2019.

# References:
- Nabighian, Misac N. (1984). “Toward a three-dimensional automatic interpretation of potential field data via generalized Hilbert transforms: Fundamental relations”. In: GEOPHYSICS 49 (6), pp. 780–786. doi: https://doi.org/10.1190/1.1441706.
- Potts, Daniel, Gabriele Steidl, and Arthur Nieslony (2004). “Fast Convolution with Radial Kernels at Nonequispaced Knots”. In: Numerische Mathematik 98, pp. 329–351. doi: https://doi.org/10.1007/s00211-004-0538-5.

# License:
Apache 2.0: https://github.com/toster90/Hilbert-transform/blob/master/LICENSE  
