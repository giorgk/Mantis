# Mantis overview
Mantis is a high-performance server application for simulating contaminant transport from diffuse pollution sources in groundwater basins. 

A detailed discussion of the background concepts of Mantis and our Non-Point Source Assessment Toolbox (NPSAT) can be found in our paper

[Kourakos, G., F.Klein, A.Cortis, and T.Harter (2012), A groundwater nonpoint source pollution modeling framework to evaluate long-term dynamics of pollutant exceedance probabilities in wells and other discharge locations, Water Resour. Res., 48, W00L13, doi:10.1029/2011WR010813](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011WR010813)

For more information, please visit our website at [gwt.ucdavis.edu](https://gwt.ucdavis.edu/) or refer to our AGU posters:

* [Kourakos, G. , T. Harter (2021) Methods and challenges of merging heterogenous data for development of large scale nitrate simulation models. Abstract H45Y-02 presented at 2021 Fall Meeting, AGU, New Orleans LA and Online Everywhere , 13-17 Dec](https://agu2021fallmeeting-agu.ipostersessions.com/Default.aspx?s=36-8E-78-92-1B-37-6E-E0-D5-8B-B6-B8-34-FE-5F-D8)
* [Cao, Z., Kourakos. G. , B.C. Jurgens, K.E. Faulkner, C.V. Henri, T. Harter (2021) Age Validation for Nonpoint Source Assessment Tool in Central Valley, California. Abstract H15A-1037 presented at 2021 Fall Meeting, AGU, New Orleans LA and Online Everywhere, 13-17 Dec](https://agu2021fallmeeting-agu.ipostersessions.com/default.aspx?s=5B-7C-25-F8-27-08-5A-CE-31-70-2D-48-5F-21-7B-24)
* [Kourakos, G., and T., Harter (2014), Simulation of nonpoint source contamination based on adaptive mesh refinement. American Geophysical Union, Fall meeting 2014, San Francisco, California, USA.](http://subsurface.gr/wp-content/uploads/2016/08/AGU_poster_2014_red.pdf)
* [Kourakos, G., and T., Harter (2013), A Validation Framework for Non-Point Source Simulation Models: Application to the Southern California Central Valley with Spatio-Temporally Heterogenous Source Rates. American Geophysical Union, Fall meeting 2013, San Francisco, California, USA.](http://subsurface.gr/wp-content/uploads/2016/08/Giorgos-AGU_2013_NPS_poster_red.pdf).

# Outline of this repository
Here, we briefly outline the structure of the Mantis repository and describe the contents of each folder.

* **CPP** Contains the source code and executables, located in `CPP/Bin`
* **MantisSA**  contains a variant of the Mantis code designed for salt simulations, incorporating feedback mechanisms from irrigated water.
* Any additional folders within the repository are typically used for development and experimentation purposes.

# Getting the code
The code is written in C++ and relies on libraries that are available for all major operating systems.

* **Windows**: Pre-compiled executables are provided  
* **Other OS**: Users will need to compile the code from source.

Detailed instructions for obtaining and building the code can be found on the [page](page) 

# Getting Started
Our [wiki](wiki) provides comprehensive documentation, including detailed guides on using Mantis.

# Who do I talk to
The appropriate contact depends on your inquiry:

* **Usage questions**: Post in the [Discussion section](https://github.com/giorgk/Mantis/discussions) for help with Mantis.
* **Code issues**: Report problems in the [Issues](https://github.com/giorgk/Mantis/issues) tracker.

For other inquiries email us at:
* gkourakos@ucdavis.edu
* thharter@ucdavis.edu
