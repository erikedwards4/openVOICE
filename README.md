# openVOICE

openVOICE: an efficient, open-source toolkit for voice features in C

Erik Edwards (erik.edwards4@gmail.com)

================================================

openVOICE is a set of command-line tools for commonly-used low-level functions in voice analytics.
The command-line programs are written in C++ with a consistent style and interface.
The low-level functions themselves are written in C, using openBLAS for best performance.
Note that this same code could be compiled with Intel MKL for a slight performance boost on Intel processors.

The interface to each C function in BLAS-like, meaning that one specifies the input and/or output dimensions,
the matrix order as row-major or column-major, and so on. Like BLAS, the use of these functions is
meant for the developer; the end user only uses the C++ command-line tools.

In earlier versions, the low-level functions were tested with highly-optimized C++ numerical libraries,
such as Armadillo (http://arma.sourceforge.net/) and Eigen3 (http://eigen.tuxfamily.org/).
However, direct use of low-level C and openBLAS was found to be fastest in all cases,
and in some cases considerably faster, because the function is tailored to the exact input setting.

The C++ command-line programs are written in a consistent style developed for command-line tools in general.
All of these command-line tools use argtable2 (http://argtable.sourceforge.net/) for parsing
inputs and option flags. All of them allow -h (--help) as a flag to give description and usage info.

The idea for this project originated with experience with openSMILE (https://www.audeering.com/opensmile/).
This is a C++ project that is heavily used in computational paralinguistics and speech processing.
However, it is no longer open-source, and despite excellent C++ software engineering overall,
it is not super optimized for efficiency and is not intuitive to use for typical use cases.

OpenVOICE is meant to be maximally efficient, very intuitive, well-documented with help information for
simple command-line tools, and to be modular such that one can use lower-level functions in novel
combinations to easily make new features. It allows input/output in several important C++ tensor formats:
Armadillo (http://arma.sourceforge.net/), ArrayFire (https://arrayfire.com/), a minimal format
for Eigen (http://eigen.tuxfamily.org/) and NumPy (https://numpy.org/).
The later means that input/output can be directly piped to snippets of Python code that use NumPy!


## Installation
```
cd /opt
git clone https://github.com/erikedwards4/openvoice
cd /opt/openvoice
make -j4
```

Each C function can also be compiled separately; see C subdirectories for details.
To make a shared object library, from /opt/openvoice do:
```
make so
```
This could be useful, for example, if trying to use the C functions in Python or one's own application.


## Dependencies
Must have argtable2 installed (available by apt-get). For Ubuntu:
```
sudo apt-get -y install libargtable2-0
```


## Usage
See each resulting command-line tool for help (use -h or --help option).
For example:
```
/opt/openvoice/stft --help
```
Or:
```
/opt/openvoice/mfcc --help
```


## Contributing
Please open an issue first to discuss what you would like to change.


## License
[BSD 3-Clause](https://choosealicense.com/licenses/bsd-3-clause/)

