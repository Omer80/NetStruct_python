# NetStruct_python
**NetStruct** is a software package for inference and analysis of population structure from genetic data, 
based on network theory methodologies. It implements analyses detailed in the following publication:

[Greenbaum G., Templeton A.R., Bar-David S. (2016).
Inference and analysis of population structure from genetic data using network theory.](http://www.genetics.org/content/202/4/1299)

This python pipeline performs the  following procedures:
- Build genetic similarity matrix from genotype data
- Perform community detection analysis for a range of edge-removal thresholds
- Perform SAD analysis for a specific edge-removal threshold

Each procedure is run independently and creates its own output files, 
which can then be used as input for the following procedure.

### Tips and recommendations
1. The first procedure, building the similarity matrix, is usually the most time-intensive procedure 
(and the only one that depends on the number of loci; see publication). 
It is therefore recommended to run this procedure once, and use the output file, 
which contains the genetic similarity matrix, to run exploratory analyses with different thresholds etc.
2. It is possible to give the output files different names than the default names. 
To do that, when typing the command line (for any procedure), 
add the desired file name as an additional parameter at the end.

## Install
The instructions below are intended for Windows 7 operating systems, 
but are similar with other operating systems (these instructions can be found online). 
They are intended for users without any experience in python, and for computers without any previous python installations.

### Installing Anaconda
Anaconda is a platform which, besides the basic python platform, has many useful packages for science application. 

In order to install the Anaconda platform:
1. Download the Anaconda installer on \url{www.continuum.io/downloads#windows}. 
2. Run the installer. Be sure to check the ``add to PATH'' box.

After installing Anaconda, python is installed and so are most packages required by **NetStruct**.

### Installing python-igraph and deepdish packages
Besides the packages in Anaconda, **NetStruct** requires two more packages: [*deepdish*](https://github.com/uchicago-cs/deepdish)
and [*python-igraph*](http://igraph.org/python/).

1. In order to install the *python-igraph* package:
  - Download the igraph-python *.whl file from [www.igraph.org/python](www.igraph.org/python) or 
  [www.lfd.uci.edu/~gohlke/pythonlibs/#python-igraph](www.lfd.uci.edu/~gohlke/pythonlibs/#python-igraph).
  Make sure to choose the python 2.7 ("cp27") and the correct version (32-bit or 64-bit, depending on the Anaconda version installed).
  - In the command prompt, move to the folder where you have downloaded the _*.whl_ file. Type
`python -m pip install SomePackage`
with the _*.whl_ file you have downloaded instead of *SomePackage*. This should install the igraph package.
2. In order to install the *deepdish* package, in the command prompt type: `python -m pip install deepdish`

## Examples
Read the NetStruct_python_manual.pdf for examples of handling the **NetStruct**.
