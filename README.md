# noise_estimation

Code for Noise Estimation for mixture of chaotic, noise and periodic motion.

# Full code in C in folder 

      estera_c
      
# Compiling from source:

      cd estera_c
      make

# Usage:

      Edit `Estera_lim.conf`
      
      Just change beginning 6 lines according to your needs:
      
        data_filenames.csv	#File name contains of file names with data (1 column starting from first raw)
        4			#No of data files
        0			#Number of lines to be ignored
        3000			#Number of data
        3000			#Size of the windows
        500			#Shift of the windows
      
      Edit `data_filenames.csv` and put in column all filenames which you would like to calculate noise
      Each file with data should have only one column of data.
      
      run `Estera`
      
# Simplified version in python

      python3 estera.py
      
# Article:


      Krzysztof Urbanowicz, Janusz A. Hołyst, Noise-level estimation of time series using coarse-grained entropy, Article in Physical Review E, May 2003.
      
      
# Errata to the article:

      In Eq. (22) should be 0.3441717, NOT 3.441717
