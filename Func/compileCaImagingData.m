% 
% compile the ca++ raw imaging data from a list of files
% 
% version 1.0
%
% Comparison list
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function CaImagingDataSet    = compileCaImagingData(CaImagingDir, CaImagingFileList)

    CaImagingDataSet   = repmat(struct('sessionName',1, 'ROI_list', 1, ...
                                'trial', 1, 'depth', 1),length(CaImagingFileList), 1);
    for nfile = 1:length(CaImagingFileList)
        
        fname               = CaImagingFileList(nfile).name;
        load([CaImagingDir fname])
        CaImagingDataSet(nfile).sessionName = CaImagingFileList(nfile).name;
        CaImagingDataSet(nfile).ROI_list    = ROI_list;
        CaImagingDataSet(nfile).trial       = trial;
        CaImagingDataSet(nfile).depth       = depth;
        CaImagingDataSet(nfile).nimage      = nimage;
    end
