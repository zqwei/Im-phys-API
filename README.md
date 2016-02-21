README -- Datasets comparison
=============================

### Current Stage

Version 2015 - 02 - 06

-   Analysis W/O *Generated Ca++ imaging dataset* (a few questions W/ generation
    equation needs solving)

-   Figures (ZW) -- Single unit analysis

    -   Figure 1: Raw activity of all neurons sorted in different ways (improved
        sorting)

    -   Figure 2: Shape of activity distribution across trial periods
        (pre-sample, sample, delay, response)

    -   Supplementary Figure 1: Selectivity over time

    -   Figure 3: Single neuron choice probability index distribution (ROC),
        across trial periods

-   Figures (ZW) -- Dimensionality reduction based analysis for collected data

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

-   Manuscript (ZW) -- *Figure legends*

Version 2015 - 02 - 11

-   Remove (ZW) -- Figure 6

-   Figures (ZW) -- Dimensionality reduction based analysis for collected data

    -   Figure 6A: Redo dPCA using 5 components

    -   Supplementary Figure 9: dPCA trajectories

-   Figures (ZW) -- Time series based analysis for collected data

    -   Figure 6B: GPFA for collected data

    -   Supplementary Figure 10: GPFA for sessions

-   Manuscript (ZW) -- *Experimental procedure* (or *Material and Methods*) for
    data analysis

    -   Preprocessing

    -   Computation of single neuron choice probability index distribution

    -   Kicking-neurons-out analysis

    -   LDA decoder and 1st-PCA

    -   Similarity matrix

    -   dPCA

Future Version

-   Add *Generated Ca++ imaging dataset*

-   Comparing in which degree the distinction among different dataset can be
    explained from *Generated Ca++ imaging dataset*

### What is this repository for?

-   Updating progress of dataset comparison

    -   Please list briefly one's progress in **Current Stage**

-   Proposing the analysis for dataset comparison

    -   Please add(+)/remove(-) the new/old analysis in **Analysis List**
        following your initials like ZW, SD, NL, TC, KS

    -   Code will be released to individuals upon requests

-   Smoothly integrating the figures and manuscript for the dataset comparison

### Dataset list

-   Spiking neurons (NL; 99 sessions) -- *\\dataset\_prepared\_2014\_04\_14...*

-   GCaMP6s Ca++ imaging (short delay; TW; 21+4 sessions) -- *an019* ; *an022*

-   GCaMP6f Ca++ imaging (short delay; TW; 1 session) -- *an039*

-   GCaMP6f Ca++ imaging (long delay; TW; 3 sessions) -- *an039*

-   \-*Generated Ca++ imaging (GCaMP6s parameters from TW; ZW)*

-   \-*Generated Ca++ imaging (no GCaMP6f parameters yet; ZW)*

### Analysis list

-   Generating simultaneously recoding dataset for spiking neurons

-   Preprocessing of datasets (used for collected dataset or simultaneously
    recoding dataset)

    -   \+(ZW): removing the units with low trials (*\< 30*) in either trial
        type condition

    -   \+(SD): spiking dataset -- removing low firing units; removing low ROC
        units

    -   \+(ZW): spiking dataset -- maximum intersection after SD prepocessing

    -   \+(ZW): Ca++ imaging dataset -- removing low ROC units

-   *Generating Ca++ imaging from spiking neurons*

    -   \+,-(ZW): generating GCaMP6s Ca++ imaging from spiking neurons according
        to dynamics and -- **generated data seems quantitively distinct from the
        real dataset**

-   Raw activity of all neurons sorted in different ways

    -   \+(SD): Sort the neurons using correlation matrix; imagesc of each
        dataset and list one neuron in each group

    -   \+(ZW): Modified the sorting algorithm

-   Shape of activity distribution across trial periods (pre-sample, sample,
    delay, response)

    -   \+(ZW)

-   Selectivity over time

    -   \+,-(ZW): two sample z-score no sorting

    -   \+(NL, SD, ZW): sorting

-   Single neuron choice probability index distribution (ROC), across trial
    periods

    -   \+(ZW)

-   Per session decodability over time

    -   \+(ZW): LDA decoder

-   Collected population decision decodability over time

    -   \+(ZW): LDA decoder

-   Saturation of decodeability with number of neurons

    -   \+(ZW): kurtosis distribution analysis

    -   \+(ZW): kicking-neurons-out analysis

-   Decoding which trial period one is in, how fast can you tell the data that
    the trial period switched

    -   \+(ZW): single LDA decoder for 4 epochs

    -   \+(ZW): three LDA decoders for adjacent epochs

    -   \+,-(ZW): single LDA decoder with time taper

-   Measures of trial-to-trial variability, by trial period

    -   \-(ZW): *not enough trial to do this*

-   Fraction of variance captured by PCA, per session, by trial period.

    -   \+(ZW): LDA-LDA similarity matrix across time per session or collected
        data

    -   \+(ZW): PCA-PCA similarity matrix across time per session or collected
        data

    -   \+(ZW): LDA-PCA similarity matrix across time per session or collected
        data

-   dPCA -- a rotation of the PCA components to capture the variances
    exclusively for temporary dynamics or choice selection

    -   \+(ZW): dPCA code from *C.K. Machens*

-   GPFA -- Time-series analysis

    -   \+(KS, SD, ZW): GPFA code from *B. Yu*


### Analysis code

All analysis code is stored in another git and is welcome to release to single
identity upon requests.

Please email ZW (*weiz (at) janelia (dot) hhmi (dot) org*) for access to
analysis code git.