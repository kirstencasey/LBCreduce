# lbcreduce

## Installation
lbcreduce can be installed in the usual way using git clone.

Packages not included in Anaconda that are needed to run lbcreduce include:
- [ccdproc](https://ccdproc.readthedocs.io/en/latest/install.html)

## What it does
The most basic idea behind lbcreduce is to take raw images (i.e. object images, flat fields, darks, bias images, etc.) and turn them into images ready for scientific analysis.

## How to run lbcreduce
It is highly recommended that you look into all the options available to you, but the simplest way to run lbcreduce is to type the following into the terminal:
```
python reduce.py --<raw_image_directory>
```
where <raw_image_directory> is some string specifying the path to the directory containing all the raw images you would like reduced as well as any flat fields, darks, bias images, etc. to be used in the reduction.

## Explanation of lbcreduce variables and options

### Bias Images (2D bias vs. zero frame)
...

### Overscan
...

### Flat Fielding
...

### Dark Current
...

### Image Stacking
....

## Primer on Image Processing
See the following: [CCD Data Reduction Guide](https://mwcraig.github.io/ccd-as-book/00-00-Preface.html)
Thanks also go to many in the Physics and Astronomy Departments at OSU for helpful discussions.
...

Kirsten Casey: casey.395@osu.edu
