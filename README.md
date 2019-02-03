# IdealObservers
Contains code and sample stimulus set for running a standard Bayesian Ideal Observer, a Region of Interest Ideal Observer, and a Foveated Ideal Observer. 

#### Conceptual Explanation of ROI and FIO

A) A flow chart for a Region of Interest Ideal Observer. (a.1) An Ideal Observer is separately run for each small 30x30px section of a face image corresponding to a center point that is sampled every 10px. (a.2) On each simulated trial, likelihoods are found for a chosen face to male be or female. The likelihoods are themselves sums of likelihoods of individual faces for each gender. (a.3) The maximum likelihood principle is used to find performance in the gender task for each separate face section and output a performance map that shows which parts of a face are the most informative for this task.
 B) A summary of the process of the computations in the FIO for two fixations. The top panels show a fixation point that is below the eyes, which is optimal in several different face discrimination tasks, including gender identification with neutral faces. The bottom panels show a fixation that is at the tip of the nose, which is the theoretical optimal fixation point when viewing happy faces in a gender identification task. (b.1-b.3), The filtering operation for a noiseless template. (b.1), A face image is conceptually divided into bins that correspond to specific contrast sensitivity functions (CSFs) as a function of retinal eccentricity. Contrast sensitivity functions that correspond to the center of fixation preserve the higher spatial frequencies (seen as a higher contrast in red in the CSF plots), while contrast sensitivity functions that are far from the fixation position act as low-pass filters and mostly leave the low spatial frequencies (seen as a low-contrast blue in the CSF plots). (b.2), The image is transformed into the frequency domain, filtered separately by each possible CSF (here only two are shown), and then transformed back into the spatial domain, resulting in a set of differently filtered images corresponding to each bin. (b.3), Corresponding bins are then extracted from the filtered images and input into a composite image that simulates foveation. The procedures in b.1–b.3 are then repeated for each of the rest of the noiseless face images, as well as for the noisy input to the model on a particular trial. A set of response variables are then calculated, from which a set of likelihoods is found of each face given the noisy image input. (b.4), A decision of which face was shown is made by taking the maximum likelihood. Across many trials, a set of proportion correct (PC) values is found, one for each fixation point, and then combined into a heatmap. iFFT, Inverse FFT.
 
<img src="/Images/FIOandROI_flow_chart.png" height="100%" width="100%"> 


#### Fitting FIO to Human Psychophysics Data to Study Human Perception

(fill in explanation here)

<img src="/Images/HumanFIO_FreeForced.png" height="100%" width="100%">

