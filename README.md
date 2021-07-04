Updated:
* MaxVol using DEIM as the initial guess
* transformation to unbounded domain

# Todo list

Before the first release
* Reimplement piecewise linear basis with optimised performance for IRT, this is already in the TT-IRT

Gradient-enhanced FTT (could be useful for LIS as well)
* Implement piecewise Hermite interpolation basis, so we can use gradient of the density function to build the FTT
* Implement gradient enhanced spectral collocation

Local computation
* Create an interface so FTT cross can use local computation to update the function evaluation at interpolation points
