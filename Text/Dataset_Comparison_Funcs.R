htmlMD <- function(images, header, width="100%"){
  htmlMD <- paste0("<table style='width:", width, "'><tr>")
  if(length(header)>1 || nchar(header)>0){
    for (h in header){
      htmlMD <- paste0(htmlMD, "<th><center>", h,"</center></th>")
    }
    htmlMD <- paste0(htmlMD, "</tr><tr>")
  }
  for (image in images){
    htmlMD <- paste0(htmlMD, "<td><img src=", image," /></td>")
  }
  htmlMD <- paste0(htmlMD, "</tr></table>")
  cat(htmlMD)
}


verticalMD <- function(images, header, width="100%"){
  htmlMD <- paste0("<table style='width:", width, "'>")
  for (nImage in 1:length(images)){
    if(length(header)>1 || nchar(header)>0){
      htmlMD <- paste0(htmlMD, "<tr><th><center>", header[nImage],"</center></th></tr>")
    }
    htmlMD <- paste0(htmlMD, "<tr><td><img src=", images[nImage]," /></td></tr>")
  }
  htmlMD <- paste0(htmlMD, '</table>')
  cat(htmlMD)
}


tmpDatasets <- c('Shuffle_Spikes', 'Shuffle_Ca_Slow_Short_Delay_Virus', 'Shuffle_Ca_Slow_Short_Delay')
tmpDatasetsName <- c('Spike', 'AAV: GP4.12', 'Transgenic: GP4.3')


includeHtmlMD <- function(ImageDir, Datasets = tmpDatasets, DatasetsName=tmpDatasetsName, fileType = 'svg'){
  ImageDir <- paste0('../Plot/', ImageDir)
  DataFiles <- Datasets
  for (nData in 1:length(DataFiles)){
    if (fileType == "pdf"){
      tmpFileNmae <- paste0(ImageDir, DataFiles[nData], '.', fileType)
      tmpSvgFile <- paste0(ImageDir, DataFiles[nData], '.svg')
      if(file.exists(tmpFileNmae)){
        system(paste('pdf2svg', tmpFileNmae, tmpSvgFile))
        system(paste('rm', tmpFileNmae))
      }
      DataFiles[nData] <- paste0(ImageDir, DataFiles[nData], '.svg')
    } # convert pdf to svg, which gives a better results
    else{
      DataFiles[nData] <- paste0(ImageDir, DataFiles[nData], '.', fileType)
    }
      
  }
  htmlMD(DataFiles, header = DatasetsName, width = '100%')
}

tmpDatasets <- c('Shuffle_Spikes', 'Shuffle_Ca_Slow_Short_Delay_Virus', 'Shuffle_Ca_Slow_Short_Delay')
tmpDatasetsName <- c('Spike', 'AAV: GP4.12', 'Transgenic: GP4.3')


includePieMD <- function(ImageDir, Datasets = tmpDatasets, DatasetsName=tmpDatasetsName, fileType = 'svg'){
  Datasets <- c(Datasets, "Label")
  ImageDir <- paste0('../Plot/', ImageDir)
  DataFiles <- Datasets
  for (nData in 1:length(DataFiles)){
    if (fileType == "pdf"){
      tmpFileNmae <- paste0(ImageDir, DataFiles[nData], '.', fileType)
      tmpSvgFile <- paste0(ImageDir, DataFiles[nData], '.svg')
      if(file.exists(tmpFileNmae)){
        system(paste('pdf2svg', tmpFileNmae, tmpSvgFile))
        system(paste('rm', tmpFileNmae))
      }
      DataFiles[nData] <- paste0(ImageDir, DataFiles[nData], '.svg')
    } # convert pdf to svg, which gives a better results
    else{
      DataFiles[nData] <- paste0(ImageDir, DataFiles[nData], '.', fileType)
    }
    
  }
  htmlMD(DataFiles, header = DatasetsName, width = '100%')
}

includeColorbarMD <- function(ImageDir, Datasets = tmpDatasets, DatasetsName=tmpDatasetsName, fileType = 'svg'){
  Datasets <- c(Datasets, "Colorbar")
  ImageDir <- paste0('../Plot/', ImageDir)
  DataFiles <- Datasets
  for (nData in 1:length(DataFiles)){
    if (fileType == "pdf"){
      tmpFileNmae <- paste0(ImageDir, DataFiles[nData], '.', fileType)
      tmpSvgFile <- paste0(ImageDir, DataFiles[nData], '.svg')
      if(file.exists(tmpFileNmae)){
        system(paste('pdf2svg', tmpFileNmae, tmpSvgFile))
        system(paste('rm', tmpFileNmae))
      }
      DataFiles[nData] <- paste0(ImageDir, DataFiles[nData], '.svg')
    } # convert pdf to svg, which gives a better results
    else{
      DataFiles[nData] <- paste0(ImageDir, DataFiles[nData], '.', fileType)
    }
    
  }
  htmlMD(DataFiles, header = DatasetsName, width = '100%')
}

includeVerticalMD <- function(ImageDir, Datasets = tmpDatasets, DatasetsName=tmpDatasetsName, fileType = 'svg'){
  ImageDir <- paste0('../Plot/', ImageDir)
  DataFiles <- Datasets
  for (nData in 1:length(DataFiles)){
    if (fileType == "pdf"){
      tmpFileNmae <- paste0(ImageDir, DataFiles[nData], '.', fileType)
      tmpSvgFile <- paste0(ImageDir, DataFiles[nData], '.svg')
      if(file.exists(tmpFileNmae)){
        system(paste('pdf2svg', tmpFileNmae, tmpSvgFile))
        system(paste('rm', tmpFileNmae))
      }
      DataFiles[nData] <- paste0(ImageDir, DataFiles[nData], '.svg')
    } # convert pdf to svg, which gives a better results
    else{
      DataFiles[nData] <- paste0(ImageDir, DataFiles[nData], '.', fileType)
    }
    
  }
  verticalMD(DataFiles, header = DatasetsName, width = '100%')
}

