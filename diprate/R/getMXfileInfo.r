getMXfileInfo <- function(top_dir='/Volumes/quaranta/VDIPRR/HTS004', toFile=FALSE)
{
	# use PlateInfo.MBK files generated during MetaXpress image exports
	# to identify sets of file names that quantify nuclei or FUCCI

	# ALGORITHM
	# 1) start from root directory (e.g. HTS001)
	# 2) find all subdirectory names recursively
	# 3) identify all subdirs containing 'Plate'
	# 4) find 'PlateInfo.MBK' file in each 'Plate' subdir
	# 5) extract information and assemble into new data.frame
	# 6) save data.frame as CSV file in root directory

	my.file.browse <- function (root=getwd(), multiple=F) {
		# .. and list.files(root)
		x <- c( dirname(normalizePath(root)), list.files(root,full.names=T) )
		isdir <- file.info(x)$isdir
		obj <- sort(isdir,index.return=T,decreasing=T)
		isdir <- obj$x
		x <- x[obj$ix]
		lbls <- sprintf('%s%s',basename(x),ifelse(isdir,'/',''))
		lbls[1] <- sprintf('../ (%s)', basename(x[1]))

		files <- c()
		sel = -1
		while ( TRUE ) {
			sel <- menu(lbls,title=sprintf('Select file(s) (0 to quit) in folder %s:',root))
			if (sel == 0 )
				break
			if (isdir[sel]) {
				# directory, browse further
				files <- c(files, my.file.browse( x[sel], multiple ))
				break
			} else {
				# file, add to list
				files <- c(files,x[sel])
				if ( !multiple )
					break
				# remove selected file from choices
				lbls <- lbls[-sel]
				x <- x[-sel]
				isdir <- isdir[-sel]
			}
		}
		return(files)
	}

	# load dataExtractFxns.r
	# found in the dipQC github repo
	# https://github.com/QuLab-VU/dipQC.git
	if(exists('parseMBK')) { message('dataExtractFxns already loaded') } else {
		cat('\n\n')
		message('Please locate the local file <dataExtractFxns.r> in the dipQC git repo')
		def			<-	my.file.browse(root='..')
		file.loc	<-	paste0(dirname(def),'/')
		if(grep('dataExtract',def)) source(def,chdir=TRUE)
	}

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
		# if image.time not present, extract time of first image and use for entiroe plate
		temp <- d$T_POSITION
		if(!is.na(image.time)) plate.info['plate_name'] <- paste0(tail(unlist(strsplit(pn,'-')),-2),collapse="-")
		plate.info['plate_name'] <- ifelse(plate.info['plate_name']=='',NA,plate.info['plate_name'])
		d <- d[,c('OBJ_SERVER_NAME','well','SOURCE_DESCRIPTION')]
		colnames(d) <- c('file_name','well','channel')
		if(is.na(image.time))
		{ 
			d$time <- min(strptime(temp, format='%Y-%m-%d %H:%M:%S'))
		} else { d$time <- image.time }
		d$time <- as.character(d$time)
		cbind(d,as.data.frame(t(plate.info)))
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
