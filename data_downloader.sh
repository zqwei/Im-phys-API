#############################################################################
# This is used to download the precompiled and raw data used in im-phys.org
# 
# author: Ziqiang Wei
# email: weiz@janelia.hhmi.org
#
#
#############################################################################

mkdir tmp
wget -O tmp/precompiled.zip https://ndownloader.figshare.com/articles/12786296/versions/1
unzip tmp/precompiled.zip -d tmp/precompiled
mkdir nonsimBenchmarks/TempDat
mkdir simModels/TempDat
cp tmp/precompiled/*.mat nonsimBenchmarks/TempDat
cp tmp/precompiled/DataListCells.mat simModels/TempDat


# download raw of ALM-6f-TG and V1-sim-TG
wget -O tmp/raw.zip https://ndownloader.figshare.com/articles/12792587/versions/1
unzip tmp/raw.zip -d rawData

# remove temporary files
rm -rf tmp
