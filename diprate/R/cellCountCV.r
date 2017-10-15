cellCountCV	<-	function(rawCVdata, normPos=1)	{
	
	log2norm <- function(count, ids, norm_idx=1, zero=log2(0.999))
	{
		# l2
		l2 <- log2(count)
			# finds time points with no cells (0), and replace it
			# with zero (log2(0.999)), so that the data can be displayed 
			# in log scale, yet easily found.
		l2[is.infinite(l2)] <- zero
	
		norm <- numeric()
		group <- as.character(unique(ids))
		for(i in group)
		{
			d <- l2[ids == i]
			norm <- append(norm, d - d[norm_idx])
		}
	
		norm
	}

	d <- rawCVdata
	d <- d[,c('Column','Row','Cell.Nucleus')]
	# rename Cell.Nucleus to Cell.count
	colnames(d)[3]	<-	'Cell.count'

	# ensure that Row and Column are char vectors
	d$Row		<-	as.character(d$Row)
	d$Column	<-	as.character(d$Column)

	# make row index of data positions for each time point
	idx			<-	grep('Timespan',d[,1])

	# extract time of acquisition for each time point
	times		<-	as.numeric(d[,2][idx])
	time.idx	<-	seq(length(times))

	# extract the number of wells we have data for
	numWells	<-	idx[2]-idx[1]-2

	# assemble data into single structure (a)
	for(tp in time.idx)
		{	
			temp			<-	d[(idx[tp]+1):(idx[tp]+numWells),]
			temp$Time		<-	times[tp]
			ifelse(tp==time.idx[1], a <- temp, a <- rbind(a,temp))
		}

	rownames(a)	<-	seq(nrow(a))

	# Generate a column for "Well" in the format RCC
	a[nchar(a$Column)==1,'Column']	<-	paste0('0',a[nchar(a$Column)==1,'Column'])
	a$Well		<-	paste0(a$Row,a$Column)
	a$Column	<-	as.integer(a$Column)

	# calculate log2 cell #
	a			<-	a[order(a$Well),]
	a			<-	a[!is.na(a$Cell.count),]

	a$nl2		<-	log2norm(a$Cell.count, a$Well, norm_idx=normPos)
	rownames(a)	<- NULL
	
	a
}