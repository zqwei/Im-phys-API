htmlMD <- function(images, header, width="100%"){
  htmlMD <- paste("<table style='width:", width, "'><tr>", sep = '')
  if(length(header)>1 || nchar(header)>0){
    for (h in header){
      htmlMD <- paste(htmlMD, "<th><center>", h,"</center></th>", sep = '')
    }
    htmlMD <- paste(htmlMD, "</tr><tr>", sep = '')
  }
  for (image in images){
    htmlMD <- paste(htmlMD, "<td><img src=", image," /></td>", sep = '')
  }
  htmlMD <- paste(htmlMD, "</tr></table>", sep = '')
  cat(htmlMD)
}


tmpDatasets <- c('Shuffle_Spikes', 'Shuffle_Ca_Slow_Short_Delay_Virus', 'Shuffle_Ca_Slow_Short_Delay')
tmpDatasetsName <- c('Spike', 'AAV: GP4.12', 'Transgenic: GP4.3')


includeHtmlMD <- function(ImageDir, Datasets = tmpDatasets, DatasetsName=tmpDatasetsName, fileType = 'svg'){
  DataFiles <- Datasets
  for (nData in 1:length(DataFiles)){
    if (fileType == "pdf"){
      tmpFileNmae <- paste(ImageDir, DataFiles[nData], '.', fileType, sep = '')
      tmpSvgFile <- paste(ImageDir, DataFiles[nData], '.svg', sep = '')
      if(file.exists(tmpFileNmae)){
        system(paste('pdf2svg', tmpFileNmae, tmpSvgFile))
        system(paste('rm', tmpFileNmae))
      }
      DataFiles[nData] <- paste(ImageDir, DataFiles[nData], '.svg', sep = '')
    } # convert pdf to svg, which gives a better results
    else{
      DataFiles[nData] <- paste(ImageDir, DataFiles[nData], '.', fileType, sep = '')
    }
      
  }
  htmlMD(DataFiles, header = DatasetsName, width = '100%')
}
