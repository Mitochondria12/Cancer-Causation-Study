# Cancer Genomics Explorer: Unveiling DNA Regions Impacting Cancer Progression
## Project Overview
In the complex landscape of cancer research, identifying DNA regions crucial to cancer development can revolutionize treatment plans and therapies. This project leverages the power of Convolutional Neural Networks (CNNs) to interpret genomic mutation frequency data, facilitating the accurate categorization of various cancer types and isolating influential DNA fragments.

## Objectives
1. **Cancer Classification**: To classify cancer samples by type, using their genomic mutation frequency as input data.
2. **Feature Selectio**n: To systematically pinpoint DNA fragments responsible for specific cancer types.
3. **Comparative Analysis**: To include healthy tissue samples for a more robust machine learning model.
## Key Features
### Mutational Frequency Analysis
Generates mean mutation frequency distributions for various cancer cohorts.
Produces individual sample mutation frequency distribution tables.
### Protein Location Identifier
Identifies proteins in the genome based on data from the blood protein atlas.
Prioritizes proteins situated in high-mutation-frequency regions for further investigation.
### Convolutional Neural Network(In Development)
Utilizes machine learning algorithms to analyze and categorize images generated from genomic mutation frequency distributions. These categorizations distinguish between various cancer types and healthy tissue, offering a nuanced approach to cancer diagnosis. Additionally, the system employs random sampling techniques to identify specific genomic regions that demonstrate increased susceptibility to particular forms of cancer progression, thereby providing avenues for targeted research and potential therapies.
## Technical Details
### Data Complexity
Incorporates both exome and intron data, often exceeding 1 million rows of mutation data points.
Analyzes data from thousands of different samples for each cancer type.
### Methodology
1. **Binary Search Analogy**: Starts with larger genome mutation fragments and iteratively focuses on smaller fragments to isolate regions of interest.
2. **Multiple CNN Models**: Deploys neural network models trained on varying depths of fragment reads.
### Feature Pruning
Utilizes feature selection techniques to eliminate irrelevant DNA fragments.
## Business Implications
Understanding mutation frequency distribution can have profound implications for personalized medicine, early detection, and targeted therapies. This project aims to push the boundaries of what we understand about the genetic basis of cancer, offering the medical and business sectors an invaluable tool for driving innovation.


