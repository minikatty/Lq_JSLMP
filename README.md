# Lq_JSLMP
One impulsive disturbance suppression algorithm for the audio signal. 
Solve an optimization in the adaptive filtering fashion using the *Adam* subgradient descent method.


#  Intro
 There are four files in the project:
 * **codes** file  provides three algorithms implementation of impulsive disturbance suppression in both speech and record music datasets; Make sure the voice toolbox has been installed before running the code.
 * **datasets** file provides the original speech and music used in our experiments and the raw noise/disturbance data extracted in the laboratory; 
 * **results** file is a collection of the processed speech and music in .wav format;
 * **synthetic_data** file provides the simulated numerical experiment data, codes, and some results. 
 # Additional Experiment: Impulsive disturbance suppression in records

To evaluate the performance of the proposed $\ell_q$-JSLMP algorithm on high-quality vinyl record music and ensure a fair comparison with the SAR algorithm, the music set used in this experiment consists of 20 tracks from https://eti.pg.edu.pl/katedra-systemow-automatyki/ICASSP2017. Being consistent with the study [1], perceptual evaluation of audio quality (PEAQ)[3] is the metric to evaluate resulting performance. However, from the perspective of sparse signal representation, $\ell_q$-JSLMP is a lossy reconstruction algorithm, and it would not be fair to directly compare it with SAR on the publicly available datasets. Because the proposed algorithm operates on both the temporal domain for mitigating impulsive disturbance and the DCT domain for reconstructing the target signal. Any small variation in the DCT coefficients of the target signal, after being synthesized using the OLA method, will cause a slight difference in the undisturbed data. Thus, we designed the following experiments to compare the performance differences between the algorithms: First, we applied the SAR and $\ell_q$-JSLMP algorithms directly to two contaminated music datasets. Second, we combined the two algorithms using the $\ell_q$-JSLMP algorithm on the results obtained from SAR. Finally, we tried to replace the data detected as impulsive disturbance by SAR with results obtained from $\ell_q$-JSLMP. The combined approaches aim to explore whether introducing sparse representations of music signals could improve recovery accuracy when dealing with missing data in music signals.

<div align=center>
<img src="https://github.com/minikatty/Lq_JSLMP/blob/main/figures/ImpPaper.png" width="400" >
</div>
<p align="center">
<small>
Fig.1 Impulsive disturbance pattern reported in [1]
</small>
</p>

The two contaminated music datasets use two real-world impulsive disturbance samples extracted from vinyl records. One is directly taken from the study by Ciołek *et al*. [1] and partially shown in Fig.1. Recalling Fig.2, the impulsive disturbances observed in vinyl records consist of short, sudden bursts of interference with a high amplitude, rather than being a collection of random outliers with background noise. These disturbances have a brief duration and undergo a rapid increase in amplitude before quickly decaying. The SAR algorithm utilizes these characteristics to detect the disturbance.

<div align=center>
<img src="https://github.com/minikatty/Lq_JSLMP/blob/main/figures/cracknoiseto.png" width="450" >
</div>
<p align="center">
<small>
Fig.2 $S\alpha S$-modeled impulsive disturbance and pattern
</small>
</p>

We also record some digital audio interference by the commercial record player Sony Hx-500 at 48 kHz. The detailed extraction steps are illustrated in Fig.3, primarily by visual inspection and comparison with the original audio. The extracted impulsive disturbances exhibit a similar pattern as observed in Fig.1. In addition, we also obtained background noise by analyzing the silent segments of the vinyl record. It is worth noting that the primary frequency components of the background noise are below 20 Hz, which is outside the range of human audible hearing. The only difference between the two synthesized contaminated music datasets is background noise. The synthesized dataset without background noise is labeled as *Dataset I*, while the other is labeled as *Dataset II*. 

<div align=center>
<img src="https://github.com/minikatty/Lq_JSLMP/blob/main/figures/Extraction.png" width="500" height="400" >
</div>
<p align="center">
<small>
Fig.3 Impulsive disturbance extraction using a record player.
</small>
</p>

The parameters of SAR, $\ell_q$-JSLMP, and the combined algorithm are listed in Table 1, and the resulting mean PEAQs are shown in Table 2. It shows that SAR achieves the best performance when algorithms are directly applied to the music datasets. While our proposed algorithm also can effectively remove the disturbance, its mean PEAQ score is not impressive. However, during experiments, it is discovered that applying $\ell_q$-JSLMP to the SAR-processed results improves the performance by about 0.3, marked by the *combined* method. Replacing the data classified as disturbance by SAR with the results obtained from $\ell_q$-JSLMP leads to a 0.03 improvement in terms of PEAQ. These demonstrate that incorporating prior sparsity information enhances restoration performance to the music signal. Finally, we observed that the $\ell_q$-norm ($q<1$), namely the regularization term, can further improve the performance of $\ell_q$-JSLMP at the cost of extra computations, therefore it is not used in this experiment.

<p align="center">
<small>
TABLE I: Parameters settings on different datasets.
</small>
</p>
<table border="1" width="500px" cellspacing="10" align="center">
<tr>
  <th align="center"> Dataset </th>
  <th align="center"> $\ell_q$-JSLMP(p,q=1.7,1) </th>
  <th align="center"> SAR </th>
  <th align="center"> Combined(No replacement) </th>	
</tr>
<tr>
  <td rowspan="2" align="center">I</td>
  <th align ="center"> $\lambda_1=0.04$</th>
  <th align ="center"> $\mu_{\alpha}=4$</th>
  <th align ="center"> $\lambda_1=0.009$</th>
</tr>
<tr>
  <th align ="center"> $\lambda_2=0.07$</th>
  <th align ="center"> $\mu_{\beta}=4.5$</th>
  <th align ="center"> $\lambda_2=0.15$</th>
</tr>
<tr>
  <td rowspan="2" align="center">II</td>
  <th align ="center"> $\lambda_1=0.04$</th>
  <th align ="center"> $\mu_{\alpha}=3.5$</th>
  <th align ="center"> $\lambda_1=0.002$</th>
</tr>
<tr>
  <th align ="center"> $\lambda_2=0.1$</th>
  <th align ="center"> $\mu_{\beta}=4$</th>
  <th align ="center"> $\lambda_2=0.15$</th>
</tr>
</table>

\begin{table}[htbp]
\small
  \caption{Mean PEAQs($\uparrow$) of algorithms on music datasets}
\centering\label{music}
\begin{tabular}[tb]{ccccc}
\toprule
{Dataset}&\multicolumn{1}{c}{Noisy}&\multicolumn{1}{c}{$\ell_q$-JSLMP}&\multicolumn{1}{c}{SAR}&\multicolumn{1}{c}{Combined}\\
\midrule
{$\mathrm{I}$}      &-3.49&-1.84&-0.47&\textbf{-0.36}\\ \hline		      	          	                     
{$\mathrm{II}$}     &-3.03&{-2.11}&-0.84&\textbf{-0.50} \\

\bottomrule
\end{tabular}
\end{table}




# Notes
In the experiment formulas, the regularization parameters $\lambda_1$ and $\lambda_2$ correspond to $eta1$ and $eta2$ in the code. If SAR code and music data are not available on this site https://eti.pg.edu.pl/katedra-systemow-automatyki/ICASSP2017, you can also download them in the  **datasets** and **codes** files in this project.
  
# Reference
The implementation of the compared algorithm is based on the work of several researchers, which was downloaded from their respective homepages. The main corresponding works are listed as follows:

[1] Ciołek, Marcin, and Maciej Niedźwiecki. "Detection of impulsive disturbances in archive audio signals." 2017 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP). IEEE, 2017.  

[2] Wen, Fei, et al. "Robust Sparse Recovery in Impulsive Noise via $\ell_p$ - $\ell_1$ Optimization." IEEE Transactions on Signal Processing 65.1 (2016): 105-118.

[3] Kabal, Peter. "An examination and interpretation of ITU-R BS. 1387: Perceptual evaluation of audio quality." TSP Lab Technical Report, Dept. Electrical & Computer Engineering, McGill University (2002): 1-89.
