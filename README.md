### Versions

#### Future Version
-   Add *Generated Ca++ imaging dataset*
-   Comparing in which degree the distinction among different dataset can be
    explained from *Generated Ca++ imaging dataset*

#### Version 2015 - 02 - 11
-   Remove -- Figure 6
-   Figures -- Dimensionality reduction based analysis for collected data
    -   Figure 6A: Redo dPCA using 5 components
    -   Supplementary Figure 9: dPCA trajectories
-   Figures -- Time series based analysis for collected data
    -   Figure 6B: GPFA for collected data
    -   Supplementary Figure 10: GPFA for sessions
-   Manuscript -- *Experimental procedure* (or *Material and Methods*) for
    data analysis
    -   Preprocessing
    -   Computation of single neuron choice probability index distribution
    -   Kicking-neurons-out analysis
    -   LDA decoder and 1st-PCA
    -   Similarity matrix
    -   dPCA

#### Version 2015 - 02 - 06

-   Analysis W/O *Generated Ca++ imaging dataset* (a few questions W/ generation
    equation needs solving)
-   Figures -- Single unit analysis
    -   Figure 1: Raw activity of all neurons sorted in different ways (improved
        sorting)
    -   Figure 2: Shape of activity distribution across trial periods
        (pre-sample, sample, delay, response)
    -   Supplementary Figure 1: Selectivity over time
    -   Figure 3: Single neuron choice probability index distribution (ROC),
        across trial periods
-   Figures -- Dimensionality reduction based analysis for collected data
    -   Figure 4: Decodability analysis
    -   Figure 4A: Collected population decision decodability over time
    -   Supplementary Figure 2: Decision decodability over time per session
    -   Figure 4B: Saturation of decodeability with number of neurons
        (kicking-neurons out result)
    -   Supplementary Figure 3: Saturation of decodeability with number of
        neurons (kurtosis distribution analysis)
    -   Figure 4C: Decoding which trial period one is in, how fast can you tell
        the data that the trial period switched (single LDA decoder for 4
        epochs)
    -   Supplementary Figure 4: single LDA decoder for 4 epochs per session
    -   Supplementary Figure 5: three LDA decoders for adjacent epochs
    -   Figure 5: LDA-LDA, PCA-PCA, LDA-PCA similarity matrix across time
    -   Supplementary Figures 6-8: LDA-LDA, PCA-PCA, LDA-PCA similarity matrix
        across time per session
    -   Figure 6: dPCA
-   Manuscript -- *Figure legends*


### Dataset list
-   Spiking neurons (NL; 99 sessions) -- *dataset\_prepared\_2014\_04\_14*
-   ~~GCaMP6s Ca++ imaging (short delay; TW; 21+4 sessions) -- *an019* ; *an022*~~
-   GCaMP6s Ca++ imaging __AAV__
-   GCaMP6s Ca++ imaging __Transgenic__
-   ~~GCaMP6f Ca++ imaging (short delay; TW; 1 session) -- *an039*~~
-   ~~GCaMP6f Ca++ imaging (long delay; TW; 3 sessions) -- *an039*~~
-   ~~Generated Ca++ imaging (GCaMP6s parameters from TW; ZW)~~
-   ~~Generated Ca++ imaging (no GCaMP6f parameters yet; ZW)~~
-   Generated Ca++ imaging from spikes
-   Generated spikes from AAV and Transgenic datasets

### Analysis list
-   Generating simultaneously recoding dataset for spiking neurons
-   Preprocessing of datasets (used for collected dataset or simultaneously
    recoding dataset)
    -   removing the units with low trials (*< 30*) in either trial
        type condition
    -   spiking dataset -- removing low firing units; removing low ROC
        units
    -   spiking dataset -- maximum intersection after SD prepocessing
    -   Ca++ imaging dataset -- removing low ROC units

-   *Generating Ca++ imaging from spiking neurons*
    -   generating GCaMP6s Ca++ imaging from spiking neurons according
        to dynamics and -- **generated data seems quantitively distinct from the
        real dataset**

-   Raw activity of all neurons sorted in different ways (Sort the neurons using correlation matrix; 
imagesc of each dataset and list one neuron in each group)
-   ~~Shape of activity distribution across trial periods (pre-sample, sample, delay, response)~~
-   Selectivity over time (two sample z-score no sorting; sorting)
-   Single neuron choice probability index distribution (ROC), across trial periods
-   ~~Per session decodability over time (LDA decoder)~~
-   Collected population decision decodability over time (LDA decoder)
-   Saturation of decodeability with number of neurons
    -   ~~kurtosis distribution analysis~~
    -   kicking-neurons-out analysis
-   Decoding which trial period one is in, how fast can you tell the data that the trial period switched
    -   single LDA decoder for 4 epochs
    -   ~~three LDA decoders for adjacent epochs~~
    -   ~~single LDA decoder with time taper~~
-   ~~Measures of trial-to-trial variability, by trial period (*not enough trial to do this*)~~
-   Fraction of variance captured by PCA, per session, by trial period.
    -   LDA-LDA
    -   PCA-PCA
    -   LDA-PCA
-   ~~dPCA -- a rotation of the PCA components to capture the variances exclusively 
for temporary dynamics or choice selection (code from *C.K. Machens*)~~
-   ~~GPFA -- Time-series analysis (GPFA code from *B. Yu*)~~
-   Explained variance from PCA

__Please email ZW (*weiz (at) janelia (dot) hhmi (dot) org*) for access to
analysis code git.__