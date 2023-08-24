# LBCreduce

## Installation
Packages needed to run LBCreduce:
- [conda](https://www.anaconda.com)
- [ccdproc](https://ccdproc.readthedocs.io/en/latest/install.html)
- [yaml](https://yaml.org)
- [astrometry.net](http://astrometry.net/use.html)
- [Source Extractor](https://www.astromatic.net/software/sextractor/)
- [imfit](https://www.mpe.mpg.de/~erwin/code/imfit/index.html)
- [pymfit](https://github.com/johnnygreco/pymfit/tree/master)
- [ArtPop](https://artpop.readthedocs.io/en/latest/), optional

LBCreduce can be installed in the usual way using git clone. Once the repository is cloned:
```
cd LBCreduce/
python setup.py install
```

## What it does
The most basic idea behind LBCreduce is to take raw images (i.e. object images, flat fields, darks, bias images, etc.) and turn them into images ready for scientific analysis.

## How to run LBCreduce
It is highly recommended that you look into all the options available to you in the configuration file lbcreduce-config.yml, but the simplest way to run LBCreduce is to type the following into the terminal:
```
cd LBCreduce/scripts
python run_reduce.py -i <image_directory>
```
Here, <image_directory> is the directory containing all the raw images needed for processing, including flat fields, bias images, object images, etc.
To see all the command-line options available to you, type:
```
python run_reduce.py -h
```
For more control over the image processing steps you'll need to review the options in the configuration file, lbcreduce-config.yml.

## Explanation of LBCreduce variables and options

### Bias images (2D bias vs. zero frame)
...

### Overscan subtraction and trimming
...

### Flat fielding
...

## Primer on Image Processing
See the following: [CCD Data Reduction Guide](https://mwcraig.github.io/ccd-as-book/00-00-Preface.html)
Thanks also go to many in the Physics and Astronomy Departments at OSU for helpful discussions.

...

Kirsten Casey: casey.395@osu.edu
