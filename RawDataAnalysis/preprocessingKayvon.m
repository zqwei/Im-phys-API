fileDirs = {'pt203', 'pt204'};
copyTarDir = '/Users/dalin/Documents/My_Research/Dataset/ALM_Svoboda_Lab_Data/Data_In_Use/Dataset_Comparison/ImagingData/delay1f2t/';

% for nDir = fileDirs
%     localDir  = ['/Volumes/My Drive/ALM_Recording_Pole_Task_Svoboda_Lab/Kayvon_ALM_GP43/' nDir{1}];
%     fileLists = dir([localDir '/*small*.mat']);
%     for nfile = 1:length(fileLists)
%         load([localDir '/' fileLists(nfile).name])
%         if isfield(dat_small, 'roi')
%             copyfile([localDir '/' fileLists(nfile).name], [copyTarDir nDir{1} '_' fileLists(nfile).name])
%         end
%     end
% end


for nDir = {'pt203', 'pt204'}
    localDir  = ['/Volumes/druckmannlab/Ziqiang/' nDir{1} '/with_traces_061118'];
    fileLists = dir([localDir '/*small*.mat']);
    for nfile = 1:length(fileLists)
        load([localDir '/' fileLists(nfile).name])
        if isfield(dat_small, 'roi')
            copyfile([localDir '/' fileLists(nfile).name], [copyTarDir nDir{1} '_' fileLists(nfile).name])
        end
    end
end

