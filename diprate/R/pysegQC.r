# rabbitMQ/celery/py-seg processing sometimes causes problems with 
# data output file (column names not added correctly) and missing well data
# 
# Identify bad data files, correct column names and generate MXtaskArgs file
# to process missing data; can do this in place on bigdata server

fixAllBadColnames <- function(file_dir,...)
{
    if(dir.exists(file_dir))
    {
        file_paths <- dir(file_dir, full.names=TRUE)
        file_paths <- file_paths[grepl('cellcount.csv',file_paths)]
        for(f in file_paths)
            fixBadColnames(f,...)
    }

}

fixBadColnames <- function(dat,cn=c("cell_count","file_name","plate_id","well"),overwrite=FALSE)
{
    if(class(dat)=='character')
    {    
        file_path <- dat
        if(file.exists(file_path))
            temp <- read.csv(file_path, as.is=TRUE)
    } else if(class(dat)=='data.frame') temp <- dat
    
    if(colnames(temp)[1] != "cell_count")
    {
        message('Bad colnames in ',tail(unlist(strsplit(file_path,'/')),1))
        colnames(temp) <- cn
        temp <- temp[apply(temp,1,function(x) !any(is.na(x))),]
        temp <- temp[order(temp$well),]
        rownames(temp) <- NULL
        if(overwrite & class(dat)=='character') write.csv(temp, file=file_path, row.names=FALSE)
    }
    invisible(temp)
}

findMissing <- function(fi,datapath="",verbose=FALSE)
{
    # function to locate all missing cell_count data and generate MXtaskArgs file 
    # for further/repeat processing using celery/rabbitMQ/py-seg

    # <fipath> should be the path to the fileInfo file in the directory <root_dir>
    # which may or may not also contain the data in the subdirectory named <Segmentation>
    # if <datapath> is specified, Segmentation directory will not be serched for data
    # if <datapath> not specified and <Segmentation> directory does not contain any 
    # files containing "cellcounts.csv", returns NA
    
    if(class(fi)!='data.frame')
    {
        # read imageFileInfo file to get list of all images, wells and plate_ids
        fipath <- ifelse(class(fi)=='character', fi, '')
        if(file.exists(fipath))
        {
            fin <- tail(unlist(strsplit(fipath,'/')),1)
            message(paste('Reading', fin ,'as imageFileInfo file'))
            root_dir <- sub(fin,"",fipath)
            fi <- read.csv(fipath, as.is=TRUE)
        }
    }
    
    if(datapath=="" | !any(grepl('cellcount.csv',dir(datapath))))
    {
        message('No data found. Must specify <datapath> containing cellcount.csv files')
        return(NA)
    } else if(datapath !="" & dir.exists(datapath) & any(grepl('cellcount.csv',dir(datapath))))
        datafns <- dir(datapath)[grepl('cellcount.csv',dir(datapath))]
    if(verbose) message(paste('Found',length(datafns),'data files'))
    plate_ids <- unique(fi$plate_id)
    missing_pids <- setdiff(as.integer(substr(datafns,6,10)),plate_ids)

    # all well names in VDIPRR data sets (308 wells; 384-well plate without outer wells)
    fixWellName <- function(wellName)
        sapply(wellName, FUN=function(x)
            ifelse(nchar(x)==2, paste0(substr(x,1,1),0,substr(x,2,2)), x))

    wells <- sapply(2:23, function(col) sapply(LETTERS[2:15], function(row) fixWellName(paste0(row,col))))
    dimnames(wells) <- NULL

    missing_data <- data.frame()
    pb <- txtProgressBar(1,length(datafns), style=3)
    for(fn in datafns)
    {
        if(!verbose) setTxtProgressBar(pb, match(fn,datafns))
        pid <- as.integer(substr(fn,6,10))
        temp <- read.csv(file.path(datapath,fn),as.is=TRUE)
#        temp <- fixBadColnames(temp)
    # will need to change if using more than single channel
        if(nrow(temp) < 308)
        {
            missing_wells <- wells[sapply(wells, function(x) !(x %in% temp$well))]
            if(verbose) 
                message(paste(length(missing_wells),'wells without cell count data in ',
                    tail(unlist(strsplit(fn,'/')),1)))
            for(w in missing_wells)
            {
                im_fn <- fi[fi$plate_id==pid & fi$well==w,'file_name']
                missing_data <- rbind(missing_data,data.frame(file_name=im_fn,plate_id=pid,well=w))
            }
        }
    }
    close(pb)
    message(paste('Found',nrow(missing_data),'images without cell counts'))
    invisible(missing_data)
}


makeMissingTaskArgs <- function(misdat,fi,defargs=list(),root_dir='/mnt/quaranta/VDIPRR/HTS006')
{

    # if value exists in defargs should use to replace values below
    # default values
    ch2_ip <- 'None'
    gp <- '/mnt/quaranta/VDIPRR/Control_images/MXgain_new.tif'
    rp <- 'False'
    wats <- 'True'
    
    for(n in names(defargs)[names(defargs) %in% c('ch2_ip','rp','wats')])
    {
        if(n == 'ch2_ip') ch2_ip <- defargs[[n]]
        if(n == 'gp') gp <- defargs[[n]]
        if(n == 'rp') rp <- defargs[[n]]
        if(n == 'wats') wats <- defargs[[n]]
    }
    
    # generate a list of arguments for MXtasksSmoother.py
    # ordered alphabetically since this is default behavior in Python pandas.DataFrame.to_csv
        
    # should verify all file_name in misdat in fileInfo (fi)
    fns <- fi$file_name
    if(!all(misdat$file_name %in% fns))
    {
        message('Not all missing data files found in fileInfo file. Quitting')
        return(NA)
    }
    
    nip <- paste(root_dir,'/Plate',as.integer(misdat$plate_id),'/',as.character(misdat$file_name),sep='')
    pid <- as.integer(misdat$plate_id)
    sp <- paste0(root_dir,'/Segmentation')
    well <- as.character(misdat$well)
    
    out <- data.frame(
        ch2_im_path=as.character(rep(ch2_ip,length(nip))),    # channel 2 image path
        gain_path=as.character(gp),                            # gain image path
        nuc_im_path=as.character(nip),                        # nuclei image path
        plate_id=as.integer(pid),                            # plate ID
        regprops=as.character(rp),                            # converted to a python boolean (e.g. True)
        save_path=as.character(sp),                            # path of directory to save file in
        well=as.character(well),                            # well name
        ws=as.character(wats))                                # converted to a python boolean (e.g. True)
    out
}


makeTaskArgs <- function(fi,root_dir,defargs=list(),...)
{
    # default values
    rp <- 'FALSE'
    ch1_name <- 'Cy3'
    ch2_name <- 'None'
    expdir <- ''
    sp <- file.path(root_dir,'Segmentation')
    owrite <- 'TRUE'
    
    # if value exists in defargs, should use to replace values below
    for(n in names(defargs)[names(defargs) %in% c('rp','regprops','ch1_name','ch2_name','expdir','gp','sp','overwrite','owrite')])
    {
        if(n %in% c('rp','regprops')) rp <- defargs[[n]]
        if(n == 'ch1_name')  ch1_name <- defargs[[n]]
        if(n == 'ch2_name') ch2_name <- defargs[[n]]
        if(n == 'expdir') expdir <- defargs[[n]]
        if(n == 'sp') sp <- defargs[[n]]
        if(n %in% c('owrite','overwrite')) owrite <- defargs[[n]]
    }

    args <- as.list(substitute(list(...)))[-1L]
    # generate a list of arguments for MXtasksSmoother.py or MXtasksMod.py (which has two extra parameters)
    # ordered alphabetically since this is default behavior in Python pandas.DataFrame.to_csv
    
    nip <- paste(root_dir,expdir,'/Plate',as.integer(fi[fi$channel==ch1_name,'plate_id']),'/',as.character(fi[fi$channel==ch1_name,'file_name']),sep='')
    if(ch2_name != 'None')
        ch2_ip <- paste(root_dir,expdir,'/Plate',as.integer(fi[fi$channel==ch2_name,'plate_id']),'/',as.character(fi[fi$channel==ch2_name,'file_name']),sep='')
    else
        ch2_ip <- rep('None',length(nip))
    pid <- as.integer(fi[fi$channel==ch1_name,'plate_id'])
    well <- as.character(fi[fi$channel==ch1_name,'well'])
    
    out <- data.frame(
        ch2_im_path=as.character(ch2_ip),               # channel 2 image path
        nuc_im_path=as.character(nip),                  # nuclei image path
        overwrite=as.character(owrite),                 # overwrite existing data, if found
        plate_id=as.integer(pid),                       # plate ID
        regprops=as.character(rp),                      # converted to a python boolean (e.g. True)
        save_path=as.character(sp),                     # path of directory to save file in
        well=as.character(well)                         # well name
    )
    if(length(args)>=1)
    {
        message(paste(length(args),'extra args included in output'))
        for(i in seq(length(args)))
        {
            out <- cbind(out,args[i])
            colnames(out)[ncol(out)] <- names(args)[i]
        }
    }
    out
}


assemPlateData <- function(segdirpath,count_re='cellcount',pids='', toFile=FALSE)
{
    out <- data.frame()
    if(pids[1] != '')
    {
        pds <- paste0(segdirpath,paste0('/Plate',pids))
        allpds <- basename(list.dirs(segdirpath, full.names=TRUE, recursive=FALSE))
        if(length(setdiff(basename(pds),allpds)) != 0  )
        {
            message(paste('Not all plate ids passed as <pids> found in',segdirpath))
            return(invisible(NA))
        }
    } else { pds <- list.dirs(segdirpath, full.names=TRUE, recursive=FALSE) }
    
    pb <- txtProgressBar(1,length(pds), style=3)
    for(p in pds) 
    {
        files <- list.files(p, full.names=TRUE)
        files <- files[grepl(count_re,files)]
        dat <- sapply(files, function(x) read.csv(x, as.is=TRUE), simplify=FALSE)
        dat <- do.call(rbind,dat)
        rownames(dat) <- NULL
        fn <- basename(p)
        # if data acquisition is not same throughout experiment, may have different number of
        # columns in py-seg output; this will add extra columns to accruing data.frame <out>
        # to match colnames of newest data file to be added
        if(ncol(dat) > ncol(out)) 
        {
            newcn <- setdiff(colnames(dat),colnames(out))
            temp <- data.frame(matrix(ncol=length(newcn),nrow=nrow(out)))
            colnames(temp) <- newcn
            out <- cbind(out,temp)
            out <- out[,colnames(dat)]
        }
        # NOTE: saved files are not modified; may have different numbers of columns across experiment
        if(toFile) write.csv(dat, file=file.path(segdirpath,paste0(fn,'_cellcounts.csv')),row.names=FALSE)
        out <- rbind(out,dat)
        setTxtProgressBar(pb, match(p,pds))
    }
    close(pb)
    invisible(out)
}


assembleCountData <- function(filepaths)
{
    readFiles <- function(fp) {out <- read.csv(fp,as.is=TRUE); out[,!grepl('X',colnames(out))]}
    out <- sapply(filepaths, readFiles, simplify=FALSE)

    # if data acquisition is not same throughout experiment, may have different number of
    # columns in py-seg output; this will add extra columns to data.frames in list <out>
    # merge list into two data.frames with common colnames
    if(!all(sapply(out, ncol)==ncol(out[[1]])))
    {
        cn <- colnames(out[sapply(out, ncol)==max(sapply(out, ncol))][[1]])
        temp1 <- do.call(rbind, out[sapply(out, ncol) != length(cn)])
        temp2 <- do.call(rbind, out[sapply(out, ncol) == length(cn)])
        newcn <- setdiff(cn,colnames(temp1))
        temp3 <- data.frame(matrix(ncol=length(newcn),nrow=nrow(temp1)))
        colnames(temp3) <- newcn
        temp1 <- cbind(temp1,temp3)
        out <- rbind(temp1,temp2)
        out <- out[,cn]
    } else { out <- do.call(rbind,out) }
    rownames(out) <- NULL
    out
}

prepPysegArgs <- function(finfo, rootdir, expdir='', ch2_name='FITC', savefi=FALSE, momlogpath='', ...)
{
    # generic function to prepare data.frame of py-seg arguments
    # finfo == fileInfo; can be path to fileInfo file, path to or data.frame of fileInfo
    if(finfo =='new')
    {
        source('getMXfileInfo.r')
        message('generating imageFileInfo file')
        # only need to generate file once
        fi <- getMXfileInfo(rootdir)
        if(momlogpath != '' & dir.exists(momlogpath) & length(dir(momlogpath)) > 0)
            fi <- tryCatch({addPlateNames(fi,momlogpath)},error=function(cond) {message('could not add plate_names'); return(fi) })

        if(savefi)
        {
            bn <- basename(rootdir)
            # should add user input request for approval to (over)write file
            write.csv(fi, file=paste(rootdir,'/',bn,'imageFileInfo.csv',sep=''),row.names=FALSE)
        }
    } else {
        if(file.exists(normalizePath(finfo)))
        { fi <- read.csv(finfo,as.is=TRUE)
        } else { message('could not find fileInfo');return(NA)}
    }
    
    makeTaskArgs(fi,root_dir=rootdir, defargs=as.list(data.frame(ch2_name=ch2_name,expdir=expdir)),...)
}
