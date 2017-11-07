getMapInfo <- function(mapName)
{
    #' Load plate map file
    #' Function to load plate map file in either \emph{csv} or Microsoft Excel formats
    #' (\emph{.xls} or \emph{.xlsx}).
    #' \emph{deprecated——will be removed}
    #'
    #' @param mapName \emph{path} to plate map file
    #'
    #' @return data.frame of plate map information
    ifelse(grepl('.xl',mapName),
        {
            mapColNames <- as.character(gdata::read.xls(mapName, nrow=1, header=FALSE, as.is=TRUE))
        },
        mapColNames <- as.character(read.csv(mapName, nrow=1, header=FALSE, as.is=TRUE)))

    mapColNames <- gsub('\xb5','micro',mapColNames)
    mapColNames <- gsub('TR:', "", mapColNames)

    ifelse(grepl('.xl',mapName),
        map <- gdata::read.xls(mapName, head=FALSE, skip=1, as.is=TRUE),
        map <- read.csv(mapName, head=FALSE, skip=1, as.is=TRUE))
    colnames(map) <- mapColNames
    map
}

addMapInfo     <- function(dfa,path.to.map)
{
    #' Add plate map annotations to data.frame
    #'
    #' Function to load plate map file, find relevant locations by well name, and add 
    #'  annotation Information to data.frame passed as argument
    #' @param dfa data.frame of data for annotation
    #' @param path.to.map Path to annotation (plate map) file
    #' 
    
    map <- getMapInfo(path.to.map)
    # map file should be data.frame with minimal colnames:
    # c('date' or 'expt.date', 'well','cell.line',)

    cn <- colnames(map)
    if(any(grepl('[Dd]ate',cn)))
    {
        date         <- unique(map[,grep('[Dd]ate',cn)[1]])
    } else     date <- substr(path.to.map,1,10)
    date         <- gsub('-','',date)
    date         <- date[!grepl("null",date)]

    wellName     <- cn[grep('[Ww]ell',cn)]
    wellNameDFA     <- colnames(dfa)[grep('[Ww]ell',colnames(dfa))]
    if(length(wellName)>1)    wellName <- wellName[nchar(wellName)==4]
    dfa$expt.date     <- date
    dfa$cell.line <- map[match(dfa[,wellNameDFA],map[,wellName]),grep('[cC]ell.line',cn)]
    dfa$drug1     <- map[match(dfa[,wellNameDFA],map[,wellName]),grep('[dD]rug1',cn)]
    dfa$drug1.conc <- map[match(dfa[,wellNameDFA],map[,wellName]),grepl('[cC]onc1',cn) & !grepl('units',cn)]
    dfa$drug1.units <- map[match(dfa[,wellNameDFA],map[,wellName]),grepl('[cC]onc1.units',cn)]
    dfa$drug2     <- map[match(dfa[,wellNameDFA],map[,wellName]),grep('[dD]rug2',cn)]
    dfa$drug2.conc <- map[match(dfa[,wellNameDFA],map[,wellName]),grepl('[cC]onc2',cn) & !grepl('units',cn)]
    dfa$drug2.units <- map[match(dfa[,wellNameDFA],map[,wellName]),grepl('[cC]onc2.units',cn)]
    dfa
}

addDrugInfo     <- function(dfa,path.to.map)
{
    #' Pointer to addMapInfo
    #' 
    do.call(addMapInfo,args=list(dfa=dfa,path.to.map=path.to.map))
}

parseImageJdata <- function(dataFol,fileName,mapName)
{
	#' Extract cell count data from ImageJ export file
    #' 
    xl <- FALSE
    if(grepl('.xl', fileName))    xl <- TRUE
    ifelse(xl, d <- gdata::read.xls(paste0(dataFol,fileName), header=FALSE),
        d <- read.csv(paste0(dataFol,fileName), as.is=TRUE, header=FALSE))
    d         <- d[,1:2]
    d[,1]     <- as.character(d[,1])
    colnames(d) <- c('name','cell.count')
    d$acq.time     <- strptime(substr(d$name, 1, 14), "%Y%m%d%H%M%S")
    well.temp <- unlist(strsplit(d$name,'-'))
    row.temp <- LETTERS[as.numeric(substr(well.temp[grep('R',well.temp)],2,3))]
    col.temp <- substr(well.temp[grep('C',well.temp)],2,3)
    d$well     <- paste0(row.temp,col.temp)
    d$rel.time <- 0
    for(w in unique(d$well))
    {
        d[d$well==w,'rel.time'] <- difftime(d[d$well==w,'acq.time'],d[d$well==w,'acq.time'][1], units='hours')
    }
    a <- aggregate(cell.count ~ well + rel.time, data=d, FUN=sum)
    a <- a[order(a$well),]
    a$nl2     <- log2norm(a$cell.count,a$well)
    map <- getMapInfo(paste0(dataFol,mapName))
    a$cellLine <- map[match(a$well,map$Well),'Description']
    a$drug     <- 'PLX4720'
    a$conc     <- map[match(a$well,map$Well),grep('PLX',colnames(map))]
    a
}

extractImageJcounts <- function(rawData)
{
    #' Extract ImageJ Counts
    #' 
    #' 
    if(class(rawData) != 'data.frame')
    {
        message('data must be a data.frame; expecting .csv file of ImageJ data')
        return(NA)
    }
    d          <- rawData
    d         <- d[,1:2]
    colnames(d) <- c('name','cell.count')
    d$acq.time     <- strptime(substr(d$name, 1, 14), "%Y%m%d%H%M%S")
    well.temp <- unlist(strsplit(d$name,'-'))
    row.temp <- LETTERS[as.numeric(substr(well.temp[grep('R',well.temp)],2,3))]
    col.temp <- substr(well.temp[grep('C',well.temp)],2,3)
    d$well     <- paste0(row.temp,col.temp)
    d$rel.time <- 0
    for(w in unique(d$well))
    {
        d[d$well==w,'rel.time'] <- difftime(d[d$well==w,'acq.time'],d[d$well==w,'acq.time'][1], units='hours')
    }
    a <- aggregate(cell.count ~ well + rel.time, data=d, FUN=sum)
    a <- a[order(a$well),]
    a$nl2     <- log2norm(a$cell.count,a$well)
    a
}

getWellRates <- function(raw, time.range=c(70,120))
{
    #' Determine rates of proliferation of cells in each well over specified time range
    #' 
    timeName <- colnames(raw)[grep('[Tt]ime', colnames(raw))]
    wellName <- colnames(raw)[grep('[Ww]ell',colnames(raw))]
    dateName <- colnames(raw)[grep('[Dd]ate',colnames(raw))]
    if(length(wellName)>1)    wellName <- wellName[nchar(wellName)==4]
    f <- formula(paste('nl2 ~ ',timeName,' * ',wellName))
    m <- lm(f, data=raw[raw[,timeName] > time.range[1] & raw[,timeName] < time.range[2],])
    wells <- unique(raw[,wellName])
    rates <- coef(m)[grep(timeName,names(coef(m)))]
    rates <- c(rates[1],rates[-1]+rates[1])
    cl     <- unique(raw$cellLine)
    expt <- ifelse(is.null(unique(raw[,dateName])), 'unknown date',unique(raw[,dateName]))
    out     <- data.frame(Well=wells, DIP=rates, cellLine=cl, Date=expt)
    rownames(out) <- NULL
    out
}

assembleData <- function(data.loc=getwd())
{
    #' Assemble cell count data with matching plate map files
    #' 
    file.names <- dir(data.loc)
    map.names <- file.names[grepl('[Mm]ap',file.names)]
    file.names <- setdiff(file.names[grepl('csv',file.names)],map.names)
    for(fn in file.names)
    {
        map.pre <- unlist(strsplit(fn,'.csv'))
        if(is.na(map.pre))
        {
            message(paste('Could not find a matching plate map for: ',fn))
            return(NA)
        } else {
            map.name <- map.names[grepl(map.pre,map.names,fixed=TRUE)]
        }
        temp <- cellCountCV(read.csv(fn,as.is=TRUE))
        temp <- addMapInfo2(temp,map.name)
        ifelse(fn==file.names[1], out <- temp,out <- rbind(out,temp))
    }
    out$uid <- paste(out$Cell.line,out$Drug1,sep='_')
    out
}



makeDil <- function(max.conc=4000, DF=4, n=8, sig=2)
{
    #' Generate a set of serially diluted values
    #'
    #' @param max.conc numeric of maximum concentration; Default value of 4000 (nanomolar)
    #' @param DF numeric of dilution factor
    #' @param n integer of total number of dilutions (including max.conc)
    #' @param sig integer of significant digits of returned values
    #' @return \emph{n} serial dilutions of \emph{max.conc} with dilution factor \emph{DF}
    #' and significant digits \emph{sig}

    signif(DF^(seq(from=log(max.conc)/log(DF), to=log(max.conc)/log(DF)-(n-1), by=-1)),sig)
}

fixWellName <- function(wellName)
    #' Fix well name
    #'
    #' Function to convert 2-character well names to 3-character names (e.g. B1 to B01)
    #' @param wellName character of well name
    #' @return 3-character well name
    sapply(wellName, FUN=function(x)
        ifelse(nchar(x)==2, paste0(substr(x,1,1),0,substr(x,2,2)), x))

fixColNum <- function(col.num)
    #' Convert integer column number into character of length 2
    #' @param col.num integer of column number
    #' @return 2-character column name
    sapply(col.num, FUN=function(x) ifelse(nchar(x)==1, paste0(0,x),x))

removeRowsByNames <- function(dat,row.names)
{
    #' Remove rows by name
    #'
    #' Remove rows from data.frame using row names
    #' @param dat data.frame
    #' @param row.names character vector
    #' @return data.frame
    #' 
    # should make this a tryCatch function
    if(!(is.data.frame(dat) & is.character(row.names)))
    {
        message('removeRowByName expects a data.frame and character')
        return(NA)
    }
    dat[!(rownames(dat) %in% row.names),]
}

insRows <- function(data, newRows, insAfter)
{
    #' Insert rows into data.frame after a specified row name
    #'
    #' Insert rows into data.frame after a specified row name
    #' @param dat data.frame
    #' @param newRows data.frame of 1 or more rows matching colnames of data
    #' @return data.frame
    #' 
    if(!(is.data.frame(data) | is.matrix(data)))
    {
        message('Data must be either matrix or data.frame')
        return(NULL)
    }

    if(!(insAfter %in% rownames(data)))
    {
        message('insAfter must be a row name in data')
        return(NULL)
    }

    upper <- data[1:(match(insAfter, rownames(data))),]
    lower <- data[(match(insAfter, rownames(data))+1):nrow(data),]
    rbind(upper,newRows,lower)
}

parseMBK <- function(mbk_file_path)
{
    #' Parse MetaXpress MBK file
    #' 
    #' This function parses MBK files, which are generated by MetaXpress during export of image files
    #' @param mbk_file_path character or path to MBK file
    #' @return data.frame of MetaXpress image file annotations
    fn <- mbk_file_path

    # first four lines are header in different format
    h <- readLines(fn,4,warn=FALSE)
    m <- read.csv(fn,sep='\t',skip=4,header=FALSE, as.is=TRUE)

    # column names must be parsed from row 3 of header
    cn <- unlist(strsplit(h[3],','))
    cn[1] <- unlist(strsplit(cn[1],'[Columns]=\"', fixed=TRUE))[2]
    cn[length(cn)] <- unlist(strsplit(tail(cn,1),'\"', fixed=TRUE))
    colnames(m) <- cn

    # extra backslashes in several character column data to be replaced by bar
    m[,sapply(m[1,],class)=='character'] <- apply(m[,sapply(m[1,],class)=='character'], 2, function(x) gsub("\\\\","|",x))
    m
}

closestTime <- function(mytime,timevec,out="pos")
{
    #' Find closest time in vector to comparator
    #' 
    #' Function to identify the value in a time series that is closest to a comparator
    #' @param mytime Comparator time value
    #' @param timevec numeric or time series of values
    #' @param out character defining the type of output (pos = position; time = value;
    #' diff = difference between closest value and /emph{mytime})
    sapply(mytime, function(mt)
    {
        dt <- abs(difftime(timevec,mt))
        r <- switch(out,
            time = timevec[which(dt == min(dt))],
            pos = which(dt == min(dt)),
            diff = min(dt)
        )
        return(r)
    })
}

swapDrugs <- function(dat, d1='drug1', d2='drug2')
{
    #' Swap data between columns using name matching
    #' 
    #' @param dat data.frame
    #' @param d1 character found via grep in one or more columns
    #' @param d2 character found via grep in one or more columns
    #' 
    #' Function generates new temporary \code{data.frame} of columns matching \emph{d2},
    #'  copies data in columns matching \emph{d1} into columns matching \emph{d2},
    #'  then replaces columns matching \emph{d1} with columns in temporary \code{data.frame}.
    #' 
    #' Useful for swapping names of drugs in drug1 and drug2 in drug combination studies.
    #' 
    
    if(class(dat) != 'data.frame')
    {
        message(cat('swapDrugs expecting a data.frame.\n Data not modified'))
        return(dat)
    }
    if(length(grep(d1,colnames(dat))) != length(grep(d2,colnames(dat))))
    {
        message(cat('number of column names matching',d1,'and',d2,'must be the same.\n Data not modified'))
        return(dat)
    }
    nd1 <- dat[,grep(d2,colnames(dat))]
    dat[,grep(d2,colnames(dat))] <- dat[,grep(d1,colnames(dat))]
    dat[,grep(d1,colnames(dat))] <- nd1
    dat
}
