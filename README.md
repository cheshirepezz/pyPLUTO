![Poster](https://github.com/cheshirepezz/pyPLUTO/blob/master/Poster.jpg "Poster")
</p>

<p align="center">
    <img alt="Size" src="https://img.shields.io/github/repo-size/cheshirepezz/pyPLUTO">
  </a>
  <img alt="Forks" src="https://img.shields.io/github/forks/cheshirepezz/pyPLUTO">
  </a>
  <img alt="Stars" src="https://img.shields.io/github/stars/cheshirepezz/pyPLUTO">
  </a>
  <img alt="Languages" src="https://img.shields.io/github/languages/count/cheshirepezz/pyPLUTO">
  </a>
  <a href="https://github.com/cheshirepezz/pyPLUTO/graphs/contributors">
    <img alt="Contributors" src="https://img.shields.io/github/contributors/cheshirepezz/pyPLUTO">
  </a>
  <img alt="MIT Licence" src="https://img.shields.io/github/license/cheshirepezz/pyPLUTO">
  </a>
  
</p>

# pyPLUTO

##### TASK: 
      Quick Tool for Visualization of PLUTO code [Mignone 2007] data

##### AUTHORS: 
      Bhargav Vaidya [University of Leeds/MPI Astronomy, Heidelberg]
      Denis Stepanovs [MPI Astronomy, Heidelberg] 
      
This is the more usable version of the Official [pyPLUTO](https://github.com/coolastro/pyPLUTO) by Bhargav Vaidya.
The code is completely written using the Python Language. Further the GUI is developed with the Tkinter Interface.


#### Features of the code: 

1. Completely based on Python and easy to work without need of any licenses
like in IDL. 
2. The GUI environment provides a tool for quick-check of data during the
simulations runs. 
3. The code is user friendly and allows the user to even do further plotting
of contours, velocity vectors
on the surface plot.
4. Also the code can read the saved user defined variables. 

#### Limitations:

1. Only 2D views are available for 3D data. 
2. The code up-till now only has visualization tools for static data. The AMR data can not be viewed with the same. 

## Installation

The sources are located on github: https://github.com/cheshirepezz/pyPLUTO/
Open a terminal and write the command below to clone pyPLUTO repo in your PC:

```
git clone https://github.com/cheshirepezz/pyPLUTO.git
cd 
```
The code can be installed using a standard procedure. This can be done by following steps:
1. Global Install:
  * Unzip the source code: ```unzip pyPLUTO.zip```
  * Enter into the directory: ```cd pyPLUTO```
  * Install the code in the default path: ```python setup.py install```

The best practice is to create your own PYTHONPATH and do a local install in the following way!

2. Local Install:
   * Create a directory where to store this module: ```mkdir MyPython_Modules•Unzip the source code: unzip pyPLUTO-1.0.zip```
   * Enter into the directory: ```cd pyPLUTO-1.0```
   * Install the code in the directory created: ```python setup.py install --prefix=<path toMyPython_Modules>```
   * Then append the following in your ".bashrc": ```export PYTHONPATH =$HOME/MyPython_Modules/lib/python<ver>/site-packages```
   
where <ver> is the python version which the user have used to install the package.After the successful installation,  the user can start using GUI application by appending the <path toGUI_pyPLUTO.py> into their PATH.

## Manual
Further details of the Installation and Getting Started can be found [here](https://raw.githubusercontent.com/coolastro/pyPLUTO/master/doc/latex/pyPLUTO.pdf)!
  
  
  
  
  
  
