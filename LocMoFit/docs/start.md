# Getting started
LocMoFit is a tool developed in MATLAB. LocMoFit composites a set of classes and functions. You can download them [here](https://www.embl.de/download/ries/LocMoFit/) or get it as a part of [SMAP](https://github.com/jries/SMAP/) (super-resolution microscopy analysis platform), a modular analysis platform for SMLM data.
LocMoFit comes with its graphic user interface only in SMAP.

LocMoFit can be used in the following ways:
* **In MATLAB environments**, as classes and functions that can be called.
* As **a plugin of SMAP** (check section _Installation_ for [SMAP](https://github.com/jries/SMAP/) on GitHub), which also has a fully functional [stand-alone version]((https://www.embl.de/download/ries/SMAPCompiled/)).

## Requirements
Requirements differ according to the two following scenarios:
### Running in MATLAB environments (with/without SMAP)
1. MATLAB 2022a (optional) and newer. Toolboxes: Optimization, Image processing, Curve fitting, Statistics, Machine Learning, and Global optimization. 
2. Mac or Windows.
3. SMAP (optional but highly recommended).

### Running in the stand-alone version of SMAP
1. Mac or Windows.
	
2. Stand-alone version of SMAP (can be downloaded from [here](https://www.embl.de/download/ries/SMAPCompiled/)
	:::{Note}
	The stand-alone version requires no MATLAB license but is limited in extendibility. Installation notes can be downloaded [here](https://www.embl.de/download/ries/SMAPCompiled/Installation_notes_SMAP_compiled.rtf).
	:::
3. MATLAB Runtime R2022a (no MATLAB license required;). 

## Installation
### With SMAP
You can access to LocMoFit by installing SMAP. Check section _Installation_ for [SMAP on GitHub](https://github.com/jries/SMAP/).

### Without SMAP
You can downloaded the zipped files of LocMoFit [here](https://www.embl.de/download/ries/LocMoFit/). To install the code, simply unzip (usually in less than 5 min) the file and add the path of the unzipped folder to MATLAB. 

## Using LocMoFit with GUI now (SMAP required)
We recommned you to follow along the tutorial {doc}`/tutorial/quickstart` to learn LocMoFit.