# Parse MetaExpress output files
# header of ATF file (Molecular Devices output format) contains 
# most of the needed annotation information

# use single output file from baseline reads as test file
#
#
# 20160829: Problem with some data files not containing plate.name in header
# Will assume if only 6 columns that plate.name is missing
#
parseATFfile <- function(path.to.file)
{
	fn <- path.to.file
	h <- read.delim(fn, nrows=5, as.is=TRUE)[4,1]
	h <- unlist(strsplit(h,'_'))
	if(length(h) == 6)
	{	names(h) <- c('expt.name','state','date','time','instrument','plate.id') } else {
		names(h) <- c('expt.name','state','date','time','plate.name','instrument','plate.id')
	}

	d <- read.csv(fn, skip=6, sep='\t', header=TRUE)

	d <- d[,1:7]	# only keep columns 1:7
	
	renameCN <- function(coln)
	{
		new.name <- switch(coln,
			Plate.ID = "plate.id",
			Well.Name = "well",
			Run.Settings.ID = "run.settings.id",
			MEASUREMENT.SET.ID = "measurement.id",
			Total.Cells..MultiWaveScoring. = "cell.count",
			Positive.W2..MultiWaveScoring. = "s_g2.count")
		ifelse(is.null(new.name),coln,new.name)
	}
	
	colnames(d) <- sapply(colnames(d), renameCN)
	

	# ensure only known columns are included
	# proper column names
	pcn <- c("plate.id","well","run.settings.id","measurement.id","cell.count","s_g2.count")
	d <- d[,colnames(d) %in% pcn]

	d$plate.name <- h['plate.name']
	d$image.time <- paste(h['date'],paste(substr(h['time'],1,2),
		substr(h['time'],3,4),substr(h['time'],5,6),sep=':'))
	d$image.time <- strptime(d$image.time, format='%Y%m%d %T')
	d
}

