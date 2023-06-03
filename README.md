# Lq_JSLMP
One impulsive disturbance suppression algorithm for audio signal. 
Solve an optimization in the adaptive filtering fashion using the *Adam* subgradient descent method.


#  Intro
 There are four files in the project:
 * **codes** file  provides three algorithms impimplementation of impulsive disturbance suppresion in both speech and record music datasets; Make sure the voice toolbox have been installed before running the code.
 * **datasets** file provides the original speech and music used in our experiments, and the raw noise/disturbance data extracted in laboratory; 
 * **results** file is a collection of the proceesed speech and music, .wav format;
 * **synthetic_data** file provides the simulated numerical experiment data, codes and some results. 
 
 # Notes
In the experiment formulas, the regularization parameters $\lambda_1$ and $\lambda_2$ correspond to $eta1$ and $eta2$ in the code. If SAR code and music data are not available on this site https://eti.pg.edu.pl/katedra-systemow-automatyki/ICASSP2017, you can also download them in the  **datasets** and **codes** files in this project.
  
# Reference
The implementation of the compared algorithm is based on the work of several researchers, which was downloaded from their respective homepages. The main corresponding works are listed as follows:

[1] Ciołek, Marcin, and Maciej Niedźwiecki. "Detection of impulsive disturbances in archive audio signals." 2017 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP). IEEE, 2017.  

[2] Wen, Fei, et al. "Robust Sparse Recovery in Impulsive Noise via $\ell_p$ - $\ell_1$ Optimization." IEEE Transactions on Signal Processing 65.1 (2016): 105-118.
