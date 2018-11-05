
<h1>Readme</h1>

<h2>Documentation</h2>
The SARCASTIC user manual can be found in the doc folder. It will tell
you how to install the suite of tools and provide a basic tutorial that 
shows you how to use SARCASTIC to simulate a SAR collection and then view the resulting image.

<h2>Dependencies</h2>
To install the SARCASTIC tools you will need the following installed and working on your system:
* cmake
* git
* gdal
* opencl
* expat
* fftw
* boost
* CGAL
* readline
* boost

<h2>Suggested Additional Packages</h2>
SARCASTIC uses a pretty basic CAD model file format called .ply (Stanford Triangle Format). You can find out more about it at http://paulbourke.net/dataformats/ply/. There are many ways to make a .ply file from a CAD model. My preferred way at the moment is to use Sketchup Make (Its free) and export the file as a collada file before conversion to a .ply file.
A useful (and pretty awesome) tool to visualise your CAD model is meshlab which can be obtained from http://www.meshlab.net.
There are some useful python scripts in the scripts folder. To use these you will need python-2.7 installed.

<h2>Compilation</h2>
SARCASTIC is compiled using the cmake build system. To build sarcastic just type the following commands from the command line:
$ cd build
$ cmake ..
$ make
$ make install

