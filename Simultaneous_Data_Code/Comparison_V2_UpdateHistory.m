% 
% Comparison based on Simultaneously recording data
% 
% -------------------------------------------------------------------------
% version 1.0
%
% Comparison list
%
% 1.  Raw activity of all neurons sorted in different ways
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
% 3.  Selectivity over time
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
% 5.  Per session decodability over time
% 6.  Collected population decision decodability over time
% 7.  Saturation of decodeability with number of neurons
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
% 9.  Measures of trial-to-trial variability, by trial period
% 10. Fraction of variance captured by PCA, per session, by trial period.
%
% -------------------------------------------------------------------------
% version 1.1
% 
% +     Save datasets while before rerun (Simultaneously recording data). 
% +     Add analysis for "Shuffled datas code"
% - 9.  Measures of trial-to-trial variability, by trial period
% - 10. Fraction of variance captured by PCA, per session, by trial period.
% + 10. Similarity of PCA and LDA coefficient vectors as function of time
% + 11. dPCA
% + 12. GPFA
% 
% -------------------------------------------------------------------------
% version 1.2
% 
% + Add Preprocessing Data Code
% - Removes all analysis (1-4) for single neurons to "Shuffled data code"
% + Comparison between shuffled data vs non-shuffled data
% + Comparison between simultaneous recording vs non-simultaneous recording
% - Replace old version of LDA computing code with one using fitcdiscr.m 
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Generate raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_0
% 
% Comparison_V1_0_1: generate data from shuffled dataset code
%
% Comparison_V1_0_2: distribution of First Sig. P Val. Time
%
% Comparison_V1_0_2 is later deployed to "Preprocessing Data Code"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - 1.  Raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - 3.  Selectivity over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  Per session decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_5
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6
%
% Comparison_V1_6_1: The same amount of units in analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.  Saturation of decodeability with number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_7
%
% Comparison_V1_7_1: The same amount of units in analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_8
%
% Comparison_V1_8_1: The same amount of units in analysis
%
% Comparison_V1_8_2: The same amount of units in analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - 9.  Measures of trial-to-trial variability, by trial period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10->9. Fraction of variance captured by PCA, per session, by trial period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_9: not used in later
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_10: LDA-LDA, PCA-PCA, LDA-PCA
%
% Comparison_V1_10_1: LDA-LDA, PCA-PCA, LDA-PCA
%
% Comparison_V1_10_1_spike: test file using spikng data only
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. dPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. GPFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_12

%%
% -------------------------------------------------------------------------
% version 2.0
%
% Comparison list
%
% 1. Decodability of trial type over time
% 2. Decodability of epochs over time
% 3. LDA-PCA correlation as indication of network dynamics
% 4. Other analysis methods
%
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Generate raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_0  : generate data from shuffled dataset code (based on 
%                    Comparison_V1_0 && Comparison_V1_0_1) based on
%                    shuffled data set from "Shuffle code"
%                    Comparison_V2_0.m; only very high ROC units are kept
%                    to reduce the number of units in each session (this is
%                    one trick to aviod overfitting)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Decodability of trial type over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Decodability of epochs over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  LDA-PCA correlation as indication of network dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Other analysis methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_4
%