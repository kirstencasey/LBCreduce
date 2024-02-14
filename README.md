# LBCreduce

## Installation
Packages needed to run LBCreduce:
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
LBCreduce has four major components: 
- **Preliminary image reduction for Large Binocular Telescope images.** In this step, LBCreduce takes raw LBC images (flat fields, bias images, and object images) and produces useable science-quality images (albeit ones that have not yet been stacked or background subtracted -- that happens later). The configuration file for this step is *lbcreduce-config.yml*. We largely followed the [CCD Data Reduction Guide](https://mwcraig.github.io/ccd-as-book/00-00-Preface.html) during development.
- **Image registration, calibration, stacking, and background subtraction.** LBCreduce uses [astrometry.net](http://astrometry.net/use.html) to register images such that each pixel is associated with the correct location on the sky to high precision. We use PanSTARRS catalogs to calibrate image magnitudes and Source Extractor to estimate image background levels. The configuration file for this step is *register_calibration-config.yml*.
- **Diffuse light modeling.** This step is intended to meausure the structural parameters and diffuse light models of dwarf galaxies, particularly dwarf elliptical and dwarf spheroidal galaxies. To accomplish this we use [imfit](https://www.mpe.mpg.de/~erwin/code/imfit/index.html) and [pymfit](https://github.com/johnnygreco/pymfit/tree/master). You are also given the option of generating artificial galaxies using [ArtPop](https://artpop.readthedocs.io/en/latest/). The configuration file for this step is *modeling-config_imfit.yml*.
- **Surface brightness fluctuation measurements.** We can use the variation in brightness between pixels in the image (aka surface brightness fluctuations) to measure distances to galaxies. This step measures the apparent SBF magnitude in your images assuming you know the structure of the diffuse light component of the galaxy. The configuration file for this step is *modeling-config_sbf.yml*. 
    
Each of these components can be used separately or in tandem, depending on your needs. 

## Getting started
The simplest way to run each of the major components of LBCreduce are shown in this section, assuming you are in the *LBCreduce/scripts/* directory. It is *highly recommended* that you look through the relevant configuration files carefully before you start. You will need to update the image directories and file names/conventions for your situation. You can also copy the original configuration files, rename them, and move them around if you wish, as long as you tell LBCreduce (shown below). If you are using the default configuration files, you don't need to provide the config file names when you run through each of these steps. 

Sample files are provided in the *LBCreduce/lbcred_demo* directory for you to test out the diffuse light modeling and the SBF steps of the pipeline. To see an example of how to use each of these steps in tandem in a more serious way, see the *background_subtraction_tests.py* and *background_subtraction_tests_targetgal.py* scripts. 

### Preliminary image reduction
Default configuration file: *lbcreduce-config.yml*

Major function used: *lbcred.reduce*, defined in *LBCreduce/lbcred/reduce.py*.
```
python run_reduce.py --config='path/to/your/config/file/your_config_file.yml'
```
### Image registration and calibration
Default configuration file: *register_calibration-config.yml*

Major functions used: *lbcred.register_images* and *lbcred.calibrate_images*, both defined in *LBCreduce/lbcred/register_calibrate.py*.
```
python register_calibrate_images.py --config='path/to/your/config/file/your_config_file.yml'
```
### Diffuse light modeling
Default configuration file: *modeling-config_imfit.yml*

Major function used: *lbcred.modeling_imfit*, defined in *LBCreduce/lbcred/modeling_imfit.py*.

See also: *LBCreduce/lbcred_demo/*
```
python run_imfit_artpop.py --config='path/to/your/config/file/your_config_file.yml'
```
### Surface brightness fluctuations
Default configuration file: *modeling-config_sbf.yml*

Major function used: *lbcred.modeling_sbf*, defined in *LBCreduce/lbcred/modeling_sbf.py*.

See also: *LBCreduce/lbcred_demo/*
```
python run_sbf.py --config='path/to/your/config/file/your_config_file.yml'
```

## Final thoughts
A full rundown of how LBCreduce works and all the options in all the configuration files is given in my dissertation, "Surface Brightness Fluctuations of Low-Mass Galaxies in the Local Volume", available at the [OhioLink EDT Center](https://etd.ohiolink.edu). 
...

Kirsten Casey: casey.395@osu.edu
