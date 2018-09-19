% 
% obtain the location of units in ca++ imaging data
% 
% version 1.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x 1 cells (yDims number of neurons) 
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function cellLocation  = getCellLocationCaTransgenetic(CaImagingDir, fname, nUnit, AP, ML, angle)

    load([CaImagingDir fname]) 
    RoiCenterPos         = ROI_list(nUnit).centerPos;
    position_from_bregma = get_cell_position_from_bregma(RoiCenterPos, [ML, AP], angle);
    cellLocation.ML      = position_from_bregma(1);
    cellLocation.AP      = position_from_bregma(2);
    depth_str_loc        = strfind(fname, 'um');
    cellLocation.depth   = str2double(fname(depth_str_loc-3:depth_str_loc-1));