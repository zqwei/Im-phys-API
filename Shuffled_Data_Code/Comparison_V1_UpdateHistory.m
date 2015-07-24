% 
% Comparison based on shuffled unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0
%
% Comparison list (adopted from "Simultaneously data code")
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
% +     Save datasets while before rerun (Simultaneously data code). 
% +     Add analysis for "Shuffled datas code"
% - 9.  Measures of trial-to-trial variability, by trial period
% - 10. Fraction of variance captured by PCA, per session, by trial period.
% + 10. Similarity of PCA and LDA coefficient vectors as function of time
% + 11. dPCA
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
% Comparison_V1_0  : Generate raw activity of all neurons
% 
% Comparison_V1_0_1: add cellinfo fields to DataSetList, including
%                    fileName, nUnit, AP_axis, ML_axis, depth, expression,
%                    cell_type information like depth range, number of
%                    trials
%
% Comparison_V1_0_2: Gather all information of each analyzed cell; print
%                    out, number of cells
%
% Comparison_V1_0_3: Check depth dependence, ML depth dependence; AP dependence
%
% Comparison_V1_0_4: Slice cubic space for Ca++ datasets in ranges of ML,
%                    AP and depths and generate and save slicing index to 
%                    DataListShuffle.mat as 'ephysCellIndex' (also append
%                    variables like 'fileList', 'fileToAnalysis' to this
%                    mat file)
% 
% Comparison_V1_0_5: Make summary plot for sliced datasets: each row stands
%                    for a oarticular datasets; first col. : p-value of
%                    trial type as a function of time and neuronal index
%                    (sort by first largest bump onset time); second col. :
%                    Fraction of cell as function of first sig. P-val time;
%                    third col. : first sig. P-val time as a function of ML
%                    and AP locations; fourth col. first sig. P-val time as
%                    a function of depth (p-value computed using temporal
%                    filter)
%
% Comparison_V1_0_6: The same as Comparison_V1_0_5, expect that p-value
%                    computed without using temporal filter
%
% Comparison_V1_0_7: The same as Comparison_V1_0_5 second col., expect that
%                    it is cdf instead of pdf
%
% Comparison_V1_0_8: Generate Ca++ raw data (not df/f, no remove of neuropil)
%
% Comparison_V2_0  : Update Comparison_V1_0 with Comparison_V1_0_1;
%                    Add function to obtain the active units in spiking
%                    datasets (using minFiringRate for each epoches as a
%                    filter: 5 Hz in this case); mini-number of analysis
%                    file is set to 20 (otherwise, there could no valid
%                    unit for short-delay GCaMP6s transgentic datasets)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_1  : 
%
% Comparison_V1_1_1: Plot mean activity imagesc with cell class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - 3.  Selectivity over time (z-score)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  Decode for trial types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_5
% 
% Comparison_V1_5_0
%
% Comparison_V1_5_1
%
% Comparison_V1_5_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Decode for epoches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.  PCA-LDA coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_7
%
% Comparison_V1_7_1: Test for spiking dataset
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.  dPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_8
%