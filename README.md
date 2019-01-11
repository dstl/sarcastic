
# Readme

## Documentation
The SARCASTIC user manual can be found in the doc folder. It will tell
you how to install the suite of tools and provide a basic tutorial that 
shows you how to use SARCASTIC to simulate a SAR collection and then view the resulting image.

## Dependencies
To install the SARCASTIC tools you will need the following installed and working on your system:

*  cmake    (Version > 3.13.2 - on some systems called cmake3)
*  git      (Version > 2.17)
*  gdal     (Version > 2.3.2_1)
*  opencl   (Version > 1.2)
*  expat    (Version > 2.2.1)
*  fftw     (Version > 3.3.8)
*  boost    (Version > 1.68.0_1)
*  CGAL     (Version > 4.13)
*  readline (Version > 7.0.5)

	
## Suggested Additional Packages
SARCASTIC uses a pretty basic CAD model file format called .ply (Stanford Triangle Format). You can find out more about it at http://paulbourke.net/dataformats/ply/. There are many ways to make a .ply file from a CAD model. My preferred way at the moment is to use Sketchup Make (Its free) and export the file as a collada file before conversion to a .ply file.
A useful (and pretty awesome) tool to visualise your CAD model is meshlab which can be obtained from http://www.meshlab.net.
There are some useful python scripts in the scripts folder. To use these you will need python-2.7 installed.

## Compilation
SARCASTIC is compiled using the cmake build system. To build sarcastic just type the following commands from the command line:

```

$ cd build
$ cmake .. (or cmake3 ..)
$ make
$ make install

```

