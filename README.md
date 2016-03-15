# List Of Calcium Imaging Deconvolution Methods

[![Gitter](https://badges.gitter.im/zqwei/Ca-Imaging-Deconv-List.svg)](https://gitter.im/zqwei/Ca-Imaging-Deconv-List?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

Here is a collection of published calcium imaging deconvolution methods and link to their codes.

A revision of each code for __basic performance test__ is marked as revision in list.

## Main goal

As a part of my thesis work (unpublished yet), we (Ziqiang Wei and Shaul Druckmann) will try to make a comparison of exisiting calcium imaging deconvolution methods, using the [published dataset](http://crcns.org/data-sets/methods/cai-1) from Karel Svoboda's lab, where both electrical and GCaMP optical responses of a single neuron are simultaneously recorded. We laid out serveral performance indice in comparison and left a note for the performance computation, which is not clear in the original codes. This comparison project is planned to be gradually invovled. In this case, one can leave notes in __Issues tracker__, if any new comparison is missing in list or any performance should be computed towards cross comparison.

> * Dataset description: Simultaneous imaging and loose-seal cell-attached electrical recordings from neurons expressing a variety of genetically encoded calcium indicators. The data is described in: Akerboom, et al JNS 2012 and Chen, et al Nature 2013. The data provides ground truth by recording electrical and GCaMP optical responses simultaneously.
> * Tested code: all code being test is highlighted in __bold__.

More inquires should be email to: weiz AT janelia DOT hhmi DOT org


## Helmchen Model
* Model:
  * Peeling
  * Parent model: 
* Main paper:
  * http://www.hifo.uzh.ch/research/helmchen/publication/grewe2010.pdf
  * 
* Code (__Matlab__): https://github.com/HelmchenLab/CalciumSim
* Original contribution: Helmchen Lab
* Revision: 

##  SMC OOPSI
* Model:
  * Parent model: 
* Main paper:
  * 
* Note:
  * change of code:
    1. __P.k = V.spikegen.EFGinv(0.01, P, V);__ to __P.k      = log(-log(1-sum(nnorm)/V.T)/V.dt);__
* Code (__Matlab__, Python): https://github.com/jovo/smc-oopsi
* Original contribution: Josh Vogelstein
* Revision:

##  Fast OOPSI
* Model:
  * Parent model: 
* Main paper:
  * 
* Note: decaying time constant parameter $\gamma = 1 - \Delta/(1.0)$ is not updated/estimated in the code.

>   1. estimating $\gamma$ is difficult
>   2. Yaski and Friedrich (2006) showed that results are somewhat robust to minor variations in time constant  

* Code (Matlab): https://github.com/jovo/fast-oopsi
* Code (__Python__): https://github.com/liubenyuan/py-oopsi
* Original contribution: Josh Vogelstein, Benyuan Liu
* Revision: https://github.com/zqwei/py-oopsi

##  Constrained Fast OOPSI
* Model: an extension of __Fast OOPSI__
	* Method extension:
			1. strict non-negative constraint of firing rate
			2. extension of  to a general AR(p) process
			3. empirical estimation of noise prior
  * Parent model: Fast OOPSI (using [__conic programming__](http://cvxopt.org/))
* Main paper:
  * 
* Code (Matlab): https://github.com/epnev/constrained-foopsi
* Code (Python): https://github.com/epnev/constrained_foopsi_python
* Original contribution: Eftychios Pnevmatikakis, Josh Merel, Losonczy Lab
* Revision:

##  Constrained Fast OOPSI (MCMC spike inference in continuous time)
* Model: 
  * Parent model: 
* Main paper:
  * 
* Code (Matlab, _beta_): https://github.com/epnev/continuous_time_ca_sampler
* Original contribution: Eftychios Pnevmatikakis, John Merel
* Revision:

## Group LASSO initialization and spatial CNMF
* Model:
  * "Group Lasso" to detect neuronal centers and activity
  * Coordinate descent gready NMF to find neuronal activity and (non-negtive) shapes, based on group lasso initialization
  * Parent model: 
* Main paper:
  * 
* Code (Matlab, Python): https://github.com/danielso/ROI_detect
* Original contribution: Daniel So
* Revision:

## Deconvolution and demixing of calcium imaging data code
* Model: An extension of __Constrained Fast OOPSI__
	* Method extension:
		1. Dual ascent method
		2. Conic programming
		3. Nonngegative Lars
	* Parent model: Constrained Fast OOPSI (using SPGL1; CVXPY)
* Main paper:
* Code (Matlab): https://github.com/epnev/ca_source_extraction
* Code (__Python__): https://github.com/agiovann/Constrained_NMF
* Code (source extraction, Python; _alpha_): https://github.com/epnev/SOURCE_EXTRACTION_PYTHON
* Original contribution: Andrea Giovannucci and Eftychios Pnevmatikakis
* Revision:

## Sequential Image Analysis
* Model: Open source package for analysise of time-series imaging data arising from fluorescence microscopy (__Losonczy Lab__)
 * Correction of motion artifacts
 * Segmentation of imaging fields into regions of interest (ROIs)
 * Extraction of dynamic signals from ROIs
 * Parent model:
  * Constrained Fast OOPSI: https://github.com/epnev/constrained_foopsi_python
* Main paper: [Kaifosh P, Zaremba J, Danielson N, and Losonczy A. SIMA: Python software for analysis of dynamic fluorescence imaging data. Frontiers in Neuroinformatics. 2014 Aug 27; 8:77. doi: 10.3389/fninf.2014.00077.](http://journal.frontiersin.org/article/10.3389/fninf.2014.00080/full)
* Code (Python): https://github.com/losonczylab/sima
* Original contribution: Losonczy Lab

## Fast rate of innovation algorithm
* Model: 
* Main paper: [Jon Oñativia, Simon R. Schultz, and Pier Luigi Dragotti, A Finite Rate of Innovation algorithm for fast and accurate spike detection from two-photon calcium imaging. (J. Neural Eng. 10 (2013) 046017)](http://stacks.iop.org/1741-2552/10/046017)
* Code (Matlab): http://www.commsp.ee.ic.ac.uk/~jo210/src/ca_transient.zip
	* File real_data.m runs the double consistency algorithm on real data and reproduces figure 7 of the journal paper.
	* File surrogate_data.m runs the double consistency algorithm on surrogate data.
	* File fig3.m reproduces figure 3 of the journal paper.
	* File fig9.m reproduces figure 9 of the journal paper.
* Original contribution: Jon Oñativia
* Revision:
