# README #

This README documents the steps required to install sarclib

## What is this repository for? ##

sarclib is a library of tools used by the Scientific Intelligence team in Dstl.

## How do I get set up? ##

### Dependencies ##

* C compiler
* Access to a local directory for installation (default : /usr/local/dstl)
* GDAL
* readline
* expat
* fftw

### Configuration ###

```
#!c

./configure
make
make install
make docs (requires doxygen)
```


### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* @CobaltGray