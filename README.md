# Cancer-Causation-Study
This project aims to identify regions of DNA involved in cancer progression.
This project will involve creating a convoluted neural network model which takes input data in the form of mutation frequency 
across the genome images for numerous cancer samples and output which cancer type the sample belongs too.
Then identifying which fragments are responsible for which cancer types by systematically iterating numerous combinations to identify regions of interest.
We can start with larger genome mutation fragments and then generate smaller fragments, this is analogous to binary search.
This will involve multiple neural network models trained on different depths of fragment read.
When the fragments become sufficiently small enough it will be possible to look at DNA within these areas.
We obviously will want to knockout fragments which have no impact, during the feature selection, we will eliminate mutation frequency for the average cancer which
are equal to zero reducing the pixel size of each image. 
This project will also need whole mutation data for patients without cancer, so that the model can identify healthy tissue patterns.

The current program generates the mean mutation frequency distribution for cancer cohorts of interest, as well as generating an individual sample mutation frequency distribution table. 

A non functional feature at the moment is the program  identifiying the location of proteins from the blood protein atlas in the genome and there corresponding mutation frequency, it selects the proteins located in the highest mutation frequency regions for further research.

The next stages of the project will be converting each individual sample mutation frequency to a image for processing.
