getMapInfo <- function(mapName)
{
    #' Load plate map file
    #' Function to load plate map file in either \emph{csv} or Microsoft Excel formats
    #'  including xls and xlsx formats.
    #' NOTE THAT THIS FUNCTION IS DEPRECATED AND WILL BE REMOVED
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
    #'  annotation information to data.frame passed as argument
    #' @param dfa data.frame of data for annotation
    #' @param path.to.map Path to annotation (plate map) file
    #' 
    #' map file should be csv file:
    #' expecting  colnames c('date' or 'expt.date', 'well', 'cell.line', 'drug1', 'drug1.conc', 'drug1.units')
    #' @return data.frame with added columns matching column names in map file.
    map <- getMapInfo(path.to.map)
    cn <- colnames(map)

    # find/make expt.date as long as a colname for 'date' does not already exists in dfa
    if(!any(grepl('[dD]ate',colnames(dfa))))
    {
        if(any(grepl('[Dd]ate',cn))) { 
            dfa$expt.date <- unique(map[,grep('[Dd]ate',cn)[1]])
        } else if(grepl("[0-9]{2}-[0-9]{2}-[0-9]{4}", path.to.map)) {
            dfa$expt.date <- substr(path.to.map,1,10) 
        } else { 
            dfa$expt.date <- NA 
        }
    }

    wellName     <- cn[grep('[Ww]ell',cn)]
    wellNameDFA     <- colnames(dfa)[grep('[Ww]ell',colnames(dfa))]
    if(length(wellName)>1)    wellName <- wellName[nchar(wellName)==4]
    
    if(any(grep('[cC]ell.line',cn))) 
        dfa$cell.line <- as.character(map[match(dfa[,wellNameDFA],map[,wellName]),grep('[cC]ell.line',cn)])
    dfa$drug1     <- as.character(map[match(dfa[,wellNameDFA],map[,wellName]),grep('[dD]rug1$',cn)])
    dfa$drug1.conc <- as.numeric(map[match(dfa[,wellNameDFA],map[,wellName]),
        which(grepl('[cC]onc1$',cn) | grepl('[dD]rug1.conc$',cn))])
    dfa$drug1.units <- as.character(map[match(dfa[,wellNameDFA],map[,wellName]),
        which(grepl('[cC]onc1.units',cn) | grepl('[dD]rug1.units',cn))])
    if(any(grepl('[dD]rug2',cn)))
    {
        dfa$drug2     <- as.character(map[match(dfa[,wellNameDFA],map[,wellName]),grep('[dD]rug2$',cn)])
        dfa$drug2.conc <- as.numeric(map[match(dfa[,wellNameDFA],map[,wellName]),
            which(grepl('[cC]onc2$',cn) | grepl('[dD]rug2.conc$',cn))])
        dfa$drug2.units <- as.character(map[match(dfa[,wellNameDFA],map[,wellName]),
            which(grepl('[cC]onc2.units',cn) | grepl('[dD]rug2.units',cn))])
    }
    dfa
}

addDrugInfo     <- function(dfa,path.to.map)
{
    #' Pointer to addMapInfo
    #' 
    do.call(addMapInfo,args=list(dfa=dfa,path.to.map=path.to.map))
}

extractImageJcounts <- function(data_file_path)
{
	#' Extract cell count data from ImageJ export file
    #' 
    #' Function extracts cell count data from ImageJ count_nuclei macro output using
    #'  Cellavista image files as input. Only first two columns are used.
    #'  First column should be Cellavista image file name. Second column should be
    #'  cell counts.
    #'  
    #' @param data_file_path path to file exported by ImageJ count object macro
    #'  with either \code{.xls[x]} or \code{.csv} file extension.
    #'  
    #' @return data.frame with column names of `well`, `cell.count` and `time`
    containsNumeric <- function(x) !all(is.na(suppressWarnings(as.numeric(x))))
    if(grepl('[\\.]xlsx*', data_file_path))
    {
        d <- readxl::readxl(data_file_path, header=FALSE)
    } else {
        d <- read.csv(data_file_path, as.is=TRUE, header=FALSE)
    }

    # check whether first row is header (not header if any value can be coerced to numeric)
    if(!containsNumeric(d[1,])) d <- data.frame(d[-1,])
    # keep only first two columns        
    d         <- d[,1:2]
    colnames(d) <- c('file.name','cell.count')
    d$cell.count <- as.integer(d$cell.count)
    d$acq.time     <- strptime(substr(d$file.name, 1, 14), "%Y%m%d%H%M%S")
    well.temp <- unlist(strsplit(d$file.name,'-'))
    row.temp <- LETTERS[as.numeric(substr(well.temp[grep('R',well.temp)],2,3))]
    col.temp <- substr(well.temp[grep('C',well.temp)],2,3)
    d$well     <- paste0(row.temp,col.temp)
    d$time <- 0
    for(w in unique(d$well))
    {
        d[d$well==w,'time'] <- difftime(d[d$well==w,'acq.time'],d[d$well==w,'acq.time'][1], units='hours')
    }
    a <- aggregate(cell.count ~ well + time, data=d, FUN=sum)
    a <- a[order(a$well),]
    a$time <- signif(a$time,3)
    rownames(a) <- NULL
    a
}


getWellRates <- function(raw, time.range=c(70,120))
{
    #' Determine rates of proliferation of cells in each well over specified time range
    #' @param raw Raw \emph{data.frame}; expecting colnames of \code{time, well, date, cell.line, nl2}
    #'  Requires normalized log2(cell.count) \code{nl2} values as a column in \code{raw}
    #'  Uses linear model fit of  nl2 ~ time within time.range
    timeName <- colnames(raw)[grep('[Tt]ime', colnames(raw))]
    wellName <- colnames(raw)[grep('[Ww]ell',colnames(raw))]
    dateName <- colnames(raw)[grep('[Dd]ate',colnames(raw))]
    cellLineName <- colnames(raw)[grepl('[cC]ell',colnames(raw)) & grepl('[lL]ine',colnames(raw))]
    if(length(wellName)>1)    wellName <- wellName[nchar(wellName)==4]
    f <- formula(paste('nl2 ~ ',timeName,' * ',wellName))
    m <- lm(f, data=raw[raw[,timeName] > time.range[1] & raw[,timeName] < time.range[2],])
    wells <- unique(raw[,wellName])
    rates <- coef(m)[grep(timeName,names(coef(m)))]
    rates <- c(rates[1],rates[-1]+rates[1])
    cl     <- unique(raw[,cellLineName])
    expt <- ifelse(is.null(unique(raw[,dateName])), 'unknown date',unique(raw[,dateName]))
    out     <- data.frame(well=wells, DIP.rate=rates, cell.line=cl, expt.date=expt)
    rownames(out) <- NULL
    out
}

assembleData <- function(data.loc=getwd())
{
    #' Assemble cell count data with matching plate map files
    #' @param data.loc character of path to data
    #' return data.frame
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
    #'  and significant digits \emph{sig}
    #' @export
    signif(DF^(seq(from=log(max.conc)/log(DF), to=log(max.conc)/log(DF)-(n-1), by=-1)),sig)
}

fixWellName <- function(wellName)
    #' Fix well name
    #'
    #' Function to convert 2-character well names to 3-character names (e.g. B1 to B01)
    #' @param wellName character of well name
    #' @return 3-character well name
    #' @export
    sapply(wellName, FUN=function(x)
        ifelse(nchar(x)==2, paste0(substr(x,1,1),0,substr(x,2,2)), x))

fixColNum <- function(col.num)
    #' Convert integer column number into character of length 2
    #' @param col.num integer of column number
    #' @return 2-character column name
    #' @export
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
    #' @export
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
    #' @export
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
    #' @export
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
    #' return numeric of `time`, `position`, or `time difference`
    #' @export
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
    #' @export
    
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

matchNum <- function(num,comp) 
{
    #' Matching to integers to characters that may or may not have an initial '0' (e.g. '01' or '1')
    #' @param num integer
    #' @param comp character
    #' @examples
    #' matchNum(1,c('01','1','010','2','02'))
    #' @export
    
    ifelse( is.integer(num) & is.character(comp),
            grepl(paste0('^0*',num,'$'),comp), 
            NA)
}


parseCVFileName <- function(fn) 
{
    #' Parse Cellavista image file names
    #' @param fn character vector of file names
    #' @export
    sapply(fn, function(x)
    {
        ftype <- strsplit(x,'.',fixed=TRUE)[[1]][2]
        o <- strsplit(x,'.',fixed=TRUE)[[1]][1]
        o <- strsplit(o,'-')[[1]]
        tim <- strptime(o[[1]],format="%Y%m%d%H%M%S")
        fnum <- o[[2]]
        row <- LETTERS[as.integer(gsub('R','',o[[3]]))]
        col <- as.integer(gsub('C','',o[[4]]))
        well <- fixWellName(paste0(row,col))
        out <- c(x,as.character(tim),fnum,row,col,well)
        names(out) <- c('file_name','time','file_num','row','col','well')
        return(out)
    })
}

getNumChannels <- function(dir_list)
{
    #' determine the number of channels imaged in a Cellavista experiment
    #' @param dir_list list of directories and files contained within them
    #' @return integer vector of length `dir_list` or 1, if all the same value 
    #' Extract all montage images and determines how many channels are present
    #' @export
    out <- unlist(lapply(dir_list, function(x) 
    {
        a <- x[grepl("1280x1280",x)]
        length(unique(unlist(sapply(strsplit(a, '_'), function(z) z[grepl('CH',z)]))))
    }))
    ifelse(all(out==out[1]),out[1],out)
}

makeCVTaskArgs <- function(datadirs, count_chan=TRUE, verbose=TRUE, nuc_index = 1, ch2_index=0)
{
    #' Make Celery/RabbitMQ task arguments for py-seg image processing of Cellavista image files
    #' @param datadirs character vector of paths to directory(ies) containing Cellavista images
    #' @param count_chan logical whether to infer the number of channels
    #' @param verbose logical whether to output information during processing
    #' @param nuc_index integer of position of channel used to image nuclei; default is 1
    #' @param ch2_index integer of position of 2nd channel to assess positivity compared to
    #'  nuclear segmentation; default is 0 (will not assess second channel, even if it exists)
    #' @return data.frame with colnames `ch2_im_path`, `nuc_im_path`, `overwrite`
    #'  `plate_id`, `regprops`, `save_path`, and `well`

    do.call(rbind, lapply(datadirs, function(mydir)
    {
        # time series are stored in directories with integer names within the top directory
        # ts_dir = time series directories
        ts_dir <- list.dirs(mydir,full.names=FALSE, recursive=FALSE)
        # exclude empty names and names that cannot be coerced into inetegers
        ts_dir <- ts_dir[!is.na(as.integer(ts_dir[ts_dir != '']))]
        ts_dir <- ts_dir[order(as.integer(ts_dir))]

        #file_list is the list of file names in each time series directory
        file_list <- lapply(ts_dir, function(x) list.files(file.path(mydir,x)))

        # remove non-image file names
        # stitched image montages contain '1280X1280' in the name
        # example: 20160513091010_CH10_C2_R2_1280x1280.jpg
        im_file_list <- lapply(file_list, function(x)
        {
            out <- x[!grepl("1280x1280",x)]
            ftype <- tolower(sapply(out, function(z) strsplit(z,'.',fixed=TRUE)[[1]][2]))
            # keep only jpg and tiff files
            out[ftype %in% c('jpg','tiff')]
        })

        # number of time points
        n <- length(im_file_list)
        im_file_list <- lapply(im_file_list, function(x) as.data.frame(t(parseCVFileName(x))))

        f <- do.call(rbind,lapply(ts_dir, function(ts)
        {
            file_path <- sapply(im_file_list[[match(ts,ts_dir)]]$file_name, function(x) file.path(mydir,ts,x))
            cbind(data.frame(path=file_path),im_file_list[[match(ts,ts_dir)]])
        }))
        rownames(f) <- NULL
        if(verbose) message(paste('found',nrow(f),'image files total in',basename(mydir)))
        # infer the number of channels at each time point by counting unique CH values in overview images (montages)
        num_chan <- getNumChannels(file_list)
        if(verbose) message(paste('Found',num_chan,'imaging channel(s) in',basename(mydir)))
        # index position of nuclear image (assuming first); maybe should make this a passed argument
        nuc_im_path <- f$path[seq(0,nrow(f)-1,num_chan)+nuc_index]
        if(verbose) message(paste(length(nuc_im_path),'nuclear images files found in',basename(mydir)))
        if(ch2_index == 0) ch2_im_path <- 'None' # default value
        if(num_chan > 1 & ch2_index !=0) ch2_im_path <- f$path[seq(0,nrow(f)-1,num_chan)+ch2_index]
        if(verbose) message(paste(length(ch2_im_path),'unique ch2 images files found in',basename(mydir)))
        well <- f$well[seq(0,nrow(f)-1,num_chan)+nuc_index]
        data.frame( ch2_im_path=ch2_im_path,
                    nuc_im_path=nuc_im_path,
                    overwrite='TRUE',
                    plate_id=basename(dirname(nuc_im_path)),
                    regprops='TRUE',
                    save_path=file.path(dirname(dirname(nuc_im_path)),'Segmentation'),
                    well=well
                    )
    }))
}

getSegDirPaths <- function(top_dir_path='.')
{
    #' Function to find \code{Segmentation} directories within a directory.
    #' @param top_dir_path path to directory to search; default is current directory
    #'  
    #' Will search within \code{top_dir_path} for directories containing '[sS]egment'
    #'  maximum depth is two levels (top directory and any directories found in \code{top_dir_path})
    #'  
    mydirs <- list.dirs(top_dir_path,recursive=FALSE,full.names=TRUE)
    segdirs <- mydirs[grepl('[Ss]egment)',mydirs)]
    if(length(segdirs) == 0)
    {
        tempdirs <- unlist(lapply(mydirs, function(x) list.dirs(x,recursive=FALSE,full.names=TRUE)))
        mydirs <- tempdirs
        segdirs <- tempdirs[grepl('[Ss]egment',tempdirs)]
    }
    if(length(segdirs)==0)
    {
        message(paste('No Segmentation directories found in',top_dir_path))
        return(NULL)
    } else { 
        return(segdirs) 
    }
}

parseIncucyteName <- function(fn)
{
    #' Extract information from Incucyte image file names
    #' @param fn \code{character} of file name
    #' @return \code{data.frame} of \emph{file_name}, \emph{expt_id} (experiment ID), 
    #'  \emph{well}, \emph{image_pos} (image position in well), \emph{image_time}
    #'  (time image was acquired).
    
    if(!is.character(fn)) stop('parseIncucyteName expects character vectors as input')
    ss <- strsplit(fn,'.',fixed=TRUE)
    file_type <- sapply(ss, function(x) x[2])
    ss2 <- sapply(ss, function(x) strsplit(x[1],'_'))
    expt_id <- sapply(ss2, function(x) x[1])
    well <- sapply(ss2, function(x) x[2])
    img_pos <- sapply(ss2, function(x) x[3])
    img_time <- sapply(ss2, function(x) 
        as.character(strptime(paste(x[4:5],collapse=''), format='%Yy%mm%dd%Hh%Mm')))
    
    data.frame( file_name=fn,
                expt_id=expt_id,
                well=well,
                image_pos=img_pos,
                image_time=as.character(img_time))
}

makeIncTaskArgs <- function(datadirs, verbose=TRUE)
{
    #' Make Celery/RabbitMQ task arguments for py-seg image processing of Incucyte image files
    #' @param fn \code{character} of file name
    #' @return \code{data.frame} of \emph{file_name}, \emph{expt_id} (experiment ID), 
    #'  \emph{well}, \emph{image_pos} (image position in well), \emph{image_time}
    #'  (time image was acquired).
    #' 
    #' This function currently does not take other arguments to generate task arguments
    #'  such as 
    
    
    do.call(rbind, lapply(datadirs, function(mydir)
    {
        #file_list is the list of file names in each experiment directory
        file_list <- list.files(mydir)

        # remove non-image file names
        im_file_list <- unlist(lapply(file_list, function(x) 
        {
            out <- x[!grepl("1280x1280",x)]
            ftype <- tolower(sapply(out, function(z) strsplit(z,'.',fixed=TRUE)[[1]][2]))
            # keep only jpg and tiff files
            out[ftype %in% c('jpg','tiff','png')]
        }))

        im_file_df <- parseIncucyteName(im_file_list)
        full_path <- normalizePath(file.path(mydir,im_file_df$file_name))
        if(verbose) message(paste('found',length(full_path),'image files total in',mydir))
        # assuming only single channel!!
        if(verbose) message(paste(length(full_path),'nuclear images files found in',mydir))
        data.frame( ch2_im_path="None", 
                    nuc_im_path=full_path, 
                    overwrite='TRUE', 
                    plate_id=0,
                    regprops='FALSE',
                    save_path=file.path(dirname(dirname(full_path)),
                        paste(basename(dirname(full_path)),'Segmentation',sep='_')),
                    well=im_file_df$well
                    )
    }))
}
