<p align="center">
  <img src="https://raw.githubusercontent.com/PackardChan/chk2021-lengthscale-dry/main/cy.cospectrum.panel.300.png">
</p>

![GitHub top language](https://img.shields.io/github/languages/top/PackardChan/chk2021-lengthscale-dry)
![GitHub repo size](https://img.shields.io/github/repo-size/PackardChan/chk2021-lengthscale-dry)
[![GitHub license](https://img.shields.io/github/license/PackardChan/chk2021-lengthscale-dry)](LICENSE)
![GitHub last commit](https://img.shields.io/github/last-commit/PackardChan/chk2021-lengthscale-dry)
<!-- zenodo doi
-->
Data displayed in Chan, Hassanzadeh & Kuang, submitted.
<!-- search in google scholar link, dash.harvard.edu
[Search Google Scholar](

http://scholar.google.com/scholar_lookup?title=Evaluating+indices+of+blocking+anticyclones+in+terms+of+their+linear+relations+with+surface+hot+extremes&author=Chan+Hassanzadeh+Kuang&publication_year=2019&journal=Geophys.+Res.+Lett.&volume=46&pages=4904–4912
)
[Google Scholar](https://scholar.google.com/scholar?as_q=Evaluating%20indices%20of%20blocking%20anticyclones%20in%20terms%20of%20their%20linear%20relations%20with%20surface%20hot%20extremes&num=10&btnG=Search+Scholar&as_occt=any&as_sauthors=Chan,&as_ylo=2019&as_allsubj=all&hl=en&c2coff=1)
Data displayed in Chan, Hassanzadeh & Kuang, 2021, JAS
include sample plot
link for badges; crontab.sh link
-->

Citation: (To be added after acceptance)

This repository includes:
- Figures `*.pdf` in the paper
- Plotting scripts `*.ncl` and necessary netCDF data files `ensemble-wise/*.nc`
- Forcing, namelist and other files to drive the model (`fms-output/`, `modelfile/`)
- Scripts to post-process model outputs

This repository does not include:
- [Software requirements](#softreq) mentioned below
- Linear response function, the square matrix

Contact: Please find my latest email on https://orcid.org/0000-0003-1843-5566

<a name="contents"></a>
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Figure-only reproduction](#figure-only-reproduction)
  - [Requirements](#requirements)
  - [Recipe](#recipe)
- [More complete reproduction](#more-complete-reproduction)
  - [Requirements](#requirements-1)
  - [Recipe](#recipe-1)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Figure-only reproduction
[Back to contents](#contents)

### Requirements

Hardware: ![cpu](https://img.shields.io/badge/cpu-1cpu*7s-yellow)

Software: Install the following **OR** `conda create --name l1_length --file requirements.txt`
- [NCAR Command Language](http://www.ncl.ucar.edu/Download/) version 6.6.2

### Recipe
- Run `make` to plot everything
- **OR** run `make fig4` to only plot Figure 4

## More complete reproduction
[Back to contents](#contents)

### Requirements

Hardware (**per ensemble**): ![disk](https://img.shields.io/badge/disk-8.2T-blue)![cpu](https://img.shields.io/badge/cpu-640cpu*1day-yellow)

<a name="softreq"></a>
Software:
- [NCAR Command Language](http://www.ncl.ucar.edu/Download/) version 6.6.2
- Linear baroclinic instability code in MATLAB provided by Paul O'Gorman
- [MATLAB](https://www.mathworks.com/products/matlab.html) R2020a
- [GFDL dry spectral dynamical core](https://www.gfdl.noaa.gov/idealized-spectral-models-quickstart/)
- [NCO](http://nco.sourceforge.net/) version 4.7.4

### Recipe
1. There may be version difference in GFDL dry spectral dynamical core. Scripts need to be adjusted accordingly.
1. Use `fms-output/*/ensemble.sh` to submit `fms-output/*/fms-runscript.sbatch`
   - [b1ctrl](fms-output/b1ctrl): The control ensemble
   - [b2kidston](fms-output/b2kidston): The K11 ensemble
   - [b10unifZeroU](fms-output/b10unifZeroU): The 1st iteration of LRF ensemble
   - [b12unifit1](fms-output/b12unifit1): The LRF ensemble (2nd and final iteration)
   - [b13resi](fms-output/b13resi): The K11&minus;LRF ensemble
   - [b14unifm3](fms-output/b14unifm3): The &minus;1.0&times;LRF ensemble
   - [b15unifp15](fms-output/b15unifp15): The 0.5&times;LRF ensemble
   - [b16unifm15](fms-output/b16unifm15): The &minus;0.5&times;LRF ensemble
1. Use [crontab-ensemble-fms-monitor.sh](crontab-ensemble-fms-monitor.sh) to further post-process model outputs
1. Run `make` to complete the process

