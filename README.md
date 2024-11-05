## Directions

Running this code requires the following dependencies:

- RStudio \>= 2023.09.1+494 with the following packages:
  - base \>= 4.3.2
  - tidyverse \>= 2.0.0
  - arsenal \>= 3.6.3
- FreeSurfer == 7.4.1
- Singularity == 3.9.8
- miniconda3 \>= 4.9.2

Scripts with .sbatch extensions are written for submission to a SLURM
batch processing system on an HPC. It is highly recommended to conduct
this analysis on an HPC. Scripts with .sh extensions are written for
Mate Desktop and can be run with either bash or zsh. QC scripts are not necessary for replicating these results but are highly recommended if you would like to repeat the analysis in another dataset.

Run each numbered script in order. Unnumbered scripts are called by the
numbered scripts and do not need to be called manually. Wait until the
script finishes before starting the next numbered script.

The original analysis was conducted on Red Hat Enterprise Linux Server
release 7.4.

The Cam-CAN data is available at https://camcan-archive.mrc-cbu.cam.ac.uk/dataaccess/

## License

Shield: [![CC BY-SA
4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by-sa/4.0/)

This work is licensed under a [Creative Commons Attribution-ShareAlike
4.0 International
License](http://creativecommons.org/licenses/by-sa/4.0/).

[![CC BY-SA
4.0](https://licensebuttons.net/l/by-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-sa/4.0/)
