## Regional Dyssynchrony Analysis ##

This handful of MATLAB functions and classes are an implementation of the regional delay time calculation presented in the following manuscript

> Suever, J. D., Fornwalt, B. K., Neuman, L. R., Delfino, J. G., Lloyd, M. S., & Oshinski, J. N. (2013). *Method to create regional mechanical dyssynchrony maps from short-axis cine steady-state free-precession images*. Journal of Magnetic Resonance Imaging, 16(4):4. http://doi.org/10.1002/jmri.24257

The goal is for the user to supply several contraction (or strain) curves from different regions over time and the algorithm will compute a reference curve from the data using QT clustering and then the delay time between any curve and this reference curve will be computed using cross-correlation analysis.

### Getting Started ###
To learn how to use the tools provided in this repository, look at [`example.m`](example.m) for a brief demonstration on how to use the tools.
Additionally, each function and class has been fully documented and the built-in `help` function can be used to read this documentation from within MATLAB.

### Mex Files ###
This project utilizes C implementations of both QT clustering and a periodic cross-correlation algorithm to speed up processing times.
You must compile these mex files prior to running the analysis. This can be done by running the `compile_mex` command from within MATLAB.



### Attribution ###
If you use this software in your academic projects, please cite the following manuscript 

> Suever, J. D., Fornwalt, B. K., Neuman, L. R., Delfino, J. G., Lloyd, M. S., & Oshinski, J. N. (2013). *Method to create regional mechanical dyssynchrony maps from short-axis cine steady-state free-precession images*. Journal of Magnetic Resonance Imaging, 16(4):4. http://doi.org/10.1002/jmri.24257

### Issues ###
Any bugs or issues can be reported using the [bug tracker](https://github.com/suever/regionalDyssynchrony/issues).

### License ###
This software is released under the 3-Clause BSD License. Please see the [LICENSE](LICENSE) file for more information.
