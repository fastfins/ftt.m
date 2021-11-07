
Updated on 07/11/2021. 

* Load the folder 'src' and its subfolders. 
* Type 'help oned' for selecting one dimensional basis.
* Type 'help FTT' for examples of building FTT. 
    * An additional example can be found in 'test_cases/example_ftt_basis.m'
* Type 'help SIRT' for examples of building SIRT. 
    * Additional examples can be found in 'test_cases/sirt_OU'
* Type 'help DIRT' for examples of building DIRT and conditional DIRT. 
    * Examples can be found in 'test_cases/dirt_double_banana'

Important change of interfaces:
* All polynomial classes have name started with Capital letters. This does not impact most of old code, but need to be careful with checking the abstract class type `Piecewise`, `Spectral` and `Oned`---those used to start with lower case iniitals.
* The DIRT class now has a simplified interface. 
    * In the simplest case, one only needs to passin the density function, the parameter dimension and the approximation domain. This automatically override the default FTT option (in DIRT class), use a default 2nd order Lagrange polynimal basis, and use the Gaussian reference. 
    * Alternatively, one can specify the density function, the parameter dimension and the approximation polynomial basis (with a given domain). This automatically override the default FTT option (in DIRT class) and use the Gaussian reference. 
    * One can also have the full specification, with polynomial basis for each layer, reference measure, FTT options, etc. 

References: 
* For SIRT, DIRT, and IRT: 
    * Cui, Dolgov and Zahm (2021). Conditional Deep Inverse Rosenblatt Transports. arXiv preprint arXiv:2106.04170.
    * Cui and Dolgov (2021). Deep composition of tensor trains using squared inverse Rosenblatt transports. Foundation of Computational Mathematics, in press. 
    * Dolgov, Anaya-Izquierdo, Fox and Scheichl (2020). Approximation and sampling of multivariate probability distributions in the tensor train decomposition. Statistics and Computing 30(3), 603-625.
* For building the tensor train using AMEN:
    * Dolgov, Savostyanov (2014). Alternating minimal energy methods for linear systems in higher dimensions. SIAM Journal on Scientific Computing 36(5), A2248-A2271.
* For functional tensor train:
    * Gorodetsky, Karaman and Marzouk (2018). A continuous analogue of the tensor-train decomposition. Computer Methods in Applied Mechanics and Engineering 347, 59-84.
    * Bigoni, Engsig-Karup, Marzouk (2016). Spectral tensor-train decomposition. SIAM Journal on Scientific Computing 38(4), A2405-A2439.
