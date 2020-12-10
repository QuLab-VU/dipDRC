getMXfileInfo <- function(top_dir='/Volumes/quaranta/VDIPRR/HTS004', toFile=FALSE)
{
	#' Parse MetaXpress PlateInfo.MBK files
	#' 
	#' Use PlateInfo.MBK files generated during MetaXpress image exports
	#'  to identify sets of file names that quantify nuclei or FUCCI
    #' @param top_dir path to directory in which PlateInfo.MBK and images are found
    #' @param toFile logical whether to write output to file
    #' 
	#' ALGORITHM
	#' \itemize{\item{Start from root directory (e.g. HTS001)},
	#' \item{Find all subdirectory names recursively},
	#' \item{Identify all subdirs containing 'Plate'},
	#' \item{Find 'PlateInfo.MBK' file in each 'Plate' subdir},
	#' \item{Extract information and assemble into new \code{data.frame}},
	#' \item{Save data.frame as CSV file in root directory}}

	getPlateInfo <- function(MBKfilePath)
	{
		d <- tryCatch(
			{
				parseMBK(paste0(MBKfilePath,'/PlateInfo.MBK'))
			},
			error=function(cond)
			{
				message(paste('error',cond,'in:',MBKfilePath))
				return(NA)
			}
		)
		if(is.na(d)) return(NA)
		pid <- gsub('Plate','',basename(MBKfilePath))
		d$well <- fixWellName(paste0(LETTERS[d$WELL_Y],d$WELL_X))
        
        # if z-stack, capture z positions
        if(any(d$Z_INDEX > 0))
        {
            ZSTACK <- TRUE
            z_pos <- as.integer(gsub('ZStep_','',d$Z_INDEX))
        }

        # info common to all images in plate
		plate.info <- unlist(strsplit(d$DIRECTORY[1],'|',fixed=TRUE))[-1]
		plate.info <- gsub('-Use Current State--','',plate.info)
		# check whether all columns are present (expecting 5)
		pinames <- c('expt_class','expt_id','plate_name','date','plate_id')
		if(length(plate.info) < 5)
		{
			message(paste('Missing information from barcode in ',pid,'; Attempting to determine missing data type.',sep=''))
			pinames <- pinames[c(	any(grepl('DRT',plate.info)),
				any(grepl('HTS',plate.info)),
				any(grepl('[A-Z][0-9]$',plate.info) | grepl('-[A-F]$',plate.info)),
				any(grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}',plate.info)),
				any(grepl('^[0-9]{5,6}',plate.info)))]
		}
		names(plate.info) <- pinames
		plate.info <- plate.info[pinames]
		
		pn <- plate.info['plate_name']
		image.time <- as.character(strptime(paste0(unlist(strsplit(pn,'-'))[1:2],collapse=""),format='%Y%m%d%H%M%S'))
		# if image.time not present, extract time of first image and use for entire plate
		temp <- d$T_POSITION
		if(!is.na(image.time)) plate.info['plate_name'] <- paste0(tail(unlist(strsplit(pn,'-')),-2),collapse="-")
		plate.info['plate_name'] <- ifelse(plate.info['plate_name']=='',NA,plate.info['plate_name'])
		d <- d[,c('OBJ_SERVER_NAME','well','SOURCE_DESCRIPTION')]
		colnames(d) <- c('file_name','well','channel')
		if(ZSTACK) d$z_pos <- z_pos

		if(is.na(image.time))
		{ 
			d$time <- min(strptime(temp, format='%Y-%m-%d %H:%M:%S'))
		} else { d$time <- image.time }
		d$time <- as.character(d$time)
		return(cbind(d,as.data.frame(t(plate.info))))
	}
	findPlateDir <- function(dirpath)
	{
		# find directories in top directory
		dl <- tryCatch({list.dirs(dirpath,recursive=FALSE)},error=NA)
		# do not include Segmentation directory if it alredy exists
		dl <- dl[!grepl('Segmentation',dl)]
		if(length(dl) != 1 && is.na(dl[1]))
		{
			message(paste('Could not find any directories in',dirpath))
			return(NA)
		}
		
		if(!any(grepl('Plate',dl))) 
		{
			# check in subdirectories
			subdir <- sapply(dl, function(x) list.dirs(x,recursive=FALSE))
			dl <- append(dl,subdir)
		}
		
		# must have 'Plate' in directory name
		return(dl[grepl('Plate',dl)])
	}
	
	message('Finding plate directories')
	dl <- unlist(findPlateDir(top_dir))

	message('Assembling plate info')
	fileinfo <- list()
	fileinfo <- lapply(dl, function(x) tryCatch({getPlateInfo(x)},error=function(cond) {return(NA)}))
	if(any(is.na(fileinfo))) message(cat('Errors found in:',paste(dl[is.na(fileinfo)],collapse='\n')))
	fileinfo <- do.call(rbind,fileinfo[!is.na(fileinfo)])

	# NOTE: 0 byte files generate a row of NA; remove
	fileinfo <- fileinfo[!is.na(fileinfo$file_name),]
	rownames(fileinfo) <- NULL
	cat('\n')
	message('Completed assembling imageFileInfo data.frame')
	if(toFile)
	{
		fn_base <- tail(unlist(strsplit(top_dir,'/')),1)
		fn <- paste('/',fn_base,'imageFileInfo.csv',sep="")
		fp <- paste0(normalizePath(top_dir),fn)
		i = 0
		if(file.exists(fp)) message('Previous fileInfo file found!')
		while(file.exists(fp))
		{
			i = i+1
			fn <- paste('/',fn_base,'imageFileInfo_',i,'.csv',sep="")
			fp <- paste0(normalizePath(top_dir),fn)
		}
		message(paste('Writing file:',fn, ' in directory:',top_dir))
		cat('\n')
		write.csv(fileinfo, file=fp, row.names=FALSE)
	}
	invisible(fileinfo)
}
