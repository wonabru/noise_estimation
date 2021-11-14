# noise_estimation

Code for Noise Estimation for mixture of chaotic, noise and periodic motion.

# Full code in C in folder 

      estera_c
      
# Compiling from source:

      cd estera_c
      make

# Usage:

      Edit `Estera_lim.conf`
      change `!NO_FILES!` to number of filename in `datafilename.dat`
      Edit `datafilenames.dat` and put in column all filenames which you would like to calculate noise
      Each file with data should have only one column of data.
      
      run `Estera_lim`
      
# Simplified version in python

      python3 estera.py
      
# Article:


      Krzysztof Urbanowicz, Janusz A. Hołyst, Noise-level estimation of time series using coarse-grained entropy, Article in Physical Review E, May 2003.
      
      
# Errata to the article:

      In Eq. (22) should be 0.3441717, NOT 3.441717
