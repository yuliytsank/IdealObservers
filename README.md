# IdealObservers
Contains code and sample stimulus set for running a standard Bayesian Ideal Observer, a Region of Interest Ideal Observer, and a Foveated Ideal Observer set up for use with two different face discrimination tasks.  

## Explanation of What an IO, ROI, and FIO Are Used For In Human Vision Research

#### IO: 

Given knowledge of all the statistics of particular task that involves a particular level of external noise (uncertainty), an Ideal Observer is able to provide a benchmark of the maximum possible performance in that task. This is usefull for lots of things, including human perception research (which is what I use it for), because it is able to help narrow down the possbible models or mechanisms that can describe how a particular task is accomplished. (see https://www.ncbi.nlm.nih.gov/pubmed/20920517 for a thorough review of how and why ideal observers have been used in vision research). 

#### ROI:

An ROI model applied to images allows one to understand which parts contain the most task-relevant information. It makes use of an Ideal Observer by applying it to small parts of the image at a time and outputs a performance map corresponding to every possible position of the center of a small window of the images used to do the task. This is useful for classification tasks where the differences betwen the classes are concentrated in certain areas of an image. For example, in a face gender discrimination task, where the faces are spatially aligned, this model can highlight the parts of a face that hold the most discriminative (most different between the two classes) information. See the plots below for a further explanation of this example.

#### FIO:

An FIO models the foveated nature of the visual system through a series of spatial filtering operations (see explanation below for more details) based on a spatially variant contrast sensitivity function (CSF) that determines which spatial frequencies are filtered out of a stimulus, depending on a chosen simulated center of fixation to that image. It then uses an optimal decision rule to decide which class was shown on a particular simulated trial, given a specified level of added uncertainty at each pixel position. It is a very useuful tool to study eye movements because it introduces differnces in spatial resolution across the simulated visual field, which results in differences in performance for tasks like a face gender discrimination task, depending on the center of a simulated fixation position. The model outputs a performance map for every possible fixation position made to an image and is able to find a theoretical optimal intial fixation position for a particular task.    

#### ROI, FIO, and Human Psychophysics Used Together:

The figure below shows an example of what kind of insight an ROI and FIO can provide when used together: 

a) The top row shows the results of a Region of Interest Ideal Observer performance
(proportion correct) map in a gender identification task for many simulated trials where the
stimulus set was only neutral faces (left), or happy faces (right). Individual observers’ average (across
trials and blocks) initial fixation positions from a Psychophysics task where observers were allowed to freely view the stimuli while performing a face gender discrimination task condition, are overlaid in blue and the white point is the average across observers (average of blue points). The ROI map for happy faces shows that there is more discriminative information in the mouth region compared to the ROI map for neutral faces. b) The bottom row shows a performance map of an FIO on simulated trials
where the stimulus set was neutral faces (left), or happy faces (right). The FIO for happy faces
shows that there is a shift downward in the theoretical optimal point of fixation from one that is
below the eyes for neutral faces to new point that is at that the tip of the nose. In addition, there
is an increase in performance at the theoretical optimal point for happy faces vs neutral faces.
The same fixation positions as in the top panel are overlaid onto the FIO maps. Performance
between the neutral and happy maps of the ROI can be directly compared, as well as
performance between the neutral and happy maps of the FIO. However, performance cannot be
directly compared between the FIO and ROI. A noise parameter was fit to each model separately
because of large differences in efficiency between the two models compared to an ideal observer.
Using the same value for both would result in ceiling or floor effects in performance.

<img src="/Images/FIO_andROI_fixs.png" height="100%" width="100%">

The plot below shows the performance of the FIO in a vertical cross section of the image for neutral face stimuli as well as happy face stimuli. In addition, performance is shown from a human psychophysics task where observers were forced to fixate 5 different points along the midline of the face. Four of those points correspond to the forehead, eyes, nose, and mouth, while the fifth one was taken from individual preferred fixation positions found in a free-viewing (freely allowed to make eye movements) version of the task. The parameters of the FIO model were fit to to human performance using the neutral-expression face stimuli in a gender discrimination task. The same parameters were then used to run the FIO model using the happy-expresson face stimuli. The results show that humans are unable to take advantage of the extra information in the mouth region (see plots above) in happy-expression stimuli for this task and their performance does not improve as it does for the FIO model. This provides may provide some insight about the kinds of face stimuli that the human visual system is tuned to in the environment because it is presumably rare for humans to identify someone's gender while a person has a happy expression. 

<img src="/Images/HumanFIO_FreeForced.png" height="100%" width="100%">



## Conceptual Explanation of ROI and FIO Algorithms

A) A flow chart for a Region of Interest Ideal Observer. (a.1) An Ideal Observer is separately run for each small 30x30px section of a face image corresponding to a center point that is sampled every 10px. (a.2) On each simulated trial, likelihoods are found for a chosen face to male be or female. The likelihoods are themselves sums of likelihoods of individual faces for each gender. (a.3) The maximum likelihood principle is used to find performance in the gender task for each separate face section and output a performance map that shows which parts of a face are the most informative for this task.
 B) A summary of the process of the computations in the FIO for two fixations. The top panels show a fixation point that is below the eyes, which is optimal in several different face discrimination tasks, including gender identification with neutral faces. The bottom panels show a fixation that is at the tip of the nose, which is the theoretical optimal fixation point when viewing happy faces in a gender identification task. (b.1-b.3), The filtering operation for a noiseless template. (b.1), A face image is conceptually divided into bins that correspond to specific contrast sensitivity functions (CSFs) as a function of retinal eccentricity. Contrast sensitivity functions that correspond to the center of fixation preserve the higher spatial frequencies (seen as a higher contrast in red in the CSF plots), while contrast sensitivity functions that are far from the fixation position act as low-pass filters and mostly leave the low spatial frequencies (seen as a low-contrast blue in the CSF plots). (b.2), The image is transformed into the frequency domain, filtered separately by each possible CSF (here only two are shown), and then transformed back into the spatial domain, resulting in a set of differently filtered images corresponding to each bin. (b.3), Corresponding bins are then extracted from the filtered images and input into a composite image that simulates foveation. The procedures in b.1–b.3 are then repeated for each of the rest of the noiseless face images, as well as for the noisy input to the model on a particular trial. A set of response variables are then calculated, from which a set of likelihoods is found of each face given the noisy image input. (b.4), A decision of which face was shown is made by taking the maximum likelihood. Across many trials, a set of proportion correct (PC) values is found, one for each fixation point, and then combined into a heatmap. iFFT, Inverse FFT.
 
<img src="/Images/FIOandROI_flow_chart.png" height="100%" width="100%"> 

## Explanation of Code

#### IO: 

##### IO.m and runIO_gender.m

- The 'runIO_gender.m' script runs the 'IO.m' (ideal observer) function for a gender discrimination task using 80 (40 male and 40 female) 500x500px neutral-expression faces that are saved in the 'NmNf_80.mat' file. The faces are resized and spatially aligned so that the eyes are set to be 2/5 of the way down from the top of the image and the chin is 1/50 of the way from the bottom of the image. 
- The script runs the ideal observer for using several different contrast values and saves results of performance and runtime in a directory called 'SavedDir_IO_gender'. 
- The 'IO.m' function simulates many trials of a gender discrimination task. On each trial, a face draw uniformly at random (of 80 in the stimulus set), has a specified level of white noise with a specific variance added at each pixel location. This then becomes a signal which can is compared to each of the 80 original templates (noise-free images) in order to find likelihoods of the signal having come from each of the 80 templates. The likelihoods are then summed within each class (40 likelihoods for male faces, and 40 for female faces), resulting in 2 class likelihoods, the maximum of which is taken to be the chosen answer (male or female) for that trial. Only the likelhoods are calculated rather than full posterior probabilites because the function takes prior probabilities to be uniform, so they act as a scalar and don't contribute unique information for each class.  
- 'IO.m' is set up in such a way that the task that is simulated can be changed by changing the stimulus set fed into the function and/or changing the number of classes in the task. For example, the same stimulus set that is used for a gender discrimination task, can be used in a (1 of 80) face identification simply by changing the parameter 'c' in the 'runIO_gender.m' script from 2 (classes) to 80. 

#### ROI:

##### IO_ROI.m and runIO_ROI_gender.m

- The 'runIO_ROI_gender.m' script runs the 'IO_ROI.m' function for the same task and stimulus set as described for the IO above.
- The script saves the results of the performance map in a directory called 'SaveDir_ROI_gender', along with intermediate variables such as the likelihoods found after a set of trials, as well as the initial parameters used to run the function. 
- The 'IO_ROI.m' function does the same thing as the 'IO.m' function, except it uses small windows of an image at a time and runs through different positions of the window centers across the stimuli. The size of each square window is currently set to 30px and the step size to sample every 10px of the input images. These parameters can be changed in lines 87-88 of 'IO_ROI.m'.
- The 'IO_ROI.m' function uses vectorized operations to simulate many trials at once. However, the amount of RAM that is available on your computer limits the number of trials that can be run at once (currently set to 10000, but can be changed as an optional parameter fed into the function). As a result, the total number of trials needed to be run are split to into sets, where each set contains no more than than maximum number of trials specified per set. The performance results from each set of trials are then averaged. 
- Similar to the IO, the task that is being simulated can be changed by changing the stimulus set and and the paramter for the specified number of classes as needed. 

#### FIO:

##### FIO_batches.m and runFIO_batches_gender.m

- The 'runFIO_batches_gender.m' script runs the 'FIO_batches.m' function for the same task and stimulus set as described for the IO above.
- The script contains paramter settings that can be changed for the spatially variant contrast sensitivity function that is used to filter the images in order to simulate a foveated visual system.  
- The 'FIO_batches.m' function is meant to be used with large stimulus sets (larger than about 15 500x500px images). The 'batches' refers to the way the function breaks up and processes the stimulus set in parts in order to maximize the use of RAM to speed up the calculations. A subset of the fixation positions that are specified by the 'fxs' and 'fys' parameters are run at a time. In addition a subset of the stimulus set is spatially filtered at and saved at a time. After the entire stimulus set is spatially filtered for a specific set of fixation positions, the function then calculates and saves the mean and covariances matrices from the filtered templates and then deletes the filtered templates themselves in order to free up space. 
- Similar to the IO, the task that is being simulated can be changed by changing the stimulus set and and the paramter for the specified number of classes as needed. 

##### FIO_smallSet.m and runFIO_smallSet.m



