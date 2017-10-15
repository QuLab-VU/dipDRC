getMapInfo <- function(mapName)
{
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
    # dfa = data for annotation
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


# right now this is specific for PLX; must modify to make usable for any drug or treatment condition
parseImageJdata <- function(dataFol,fileName,mapName)        
{
    source('log2norm.R')
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


addDrugInfo     <- function(dfa,path.to.map)
{
    # dfa = data for annotation
    map <- getMapInfo(path.to.map)
    colnames(map) <- tolower(colnames(map))
    cn <- colnames(map)
    wellName     <- cn[grep('[Ww]ell',cn)]
    wellNameDFA     <- colnames(dfa)[grep('[Ww]ell',colnames(dfa))]
    dfa$drug1     <- map[match(dfa[,wellNameDFA],map[,wellName]),'drug1']
    dfa$drug1.conc <- map[match(dfa[,wellNameDFA],map[,wellName]),'drug1.conc']
    dfa$drug1.units <- map[match(dfa[,wellNameDFA],map[,wellName]),'drug1.units']
    dfa
}


makeDil <- function(max.conc=4000, DF=4, n=8, sig=2)
{
    # assuming max.conc in nanomolar, but any value should work.
    # max.conc = maximum concentration
    # DF = dilution factor
    # n = total number of dilutions (including max.conc)
    # sig = significant digits of output

    signif(DF^(seq(from=log(max.conc)/log(DF), to=log(max.conc)/log(DF)-(n-1), by=-1)),sig)
}

fixWellName <- function(wellName)
    sapply(wellName, FUN=function(x)
        ifelse(nchar(x)==2, paste0(substr(x,1,1),0,substr(x,2,2)), x))

fixColNum <- function(col.num)
    sapply(col.num, FUN=function(x) ifelse(nchar(x)==1, paste0(0,x),x))

removeRowsByNames <- function(dat,row.names)
{
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
# MBK files are generated by MetaXpress during export of image files
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
