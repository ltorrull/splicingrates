# inputs:
args = commandArgs(trailingOnly = T)
JUNCTIONFILE = args[1] # file with junction counts, Dprime, and Rprime
TXNRATE = as.numeric(args[2]) #txn rate in nt/min
#OUTDIR=args[3]


##### functions ######
# jointly model half-lives from all timepoints
sumsqequationsolve <- function(atts, txnrate){
	ntimes = length(atts)/2
	# first half of vector is D'
	D_prime = atts[1 : ntimes]
	# second half of vector is R'
	R_prime = atts[(ntimes+1) : length(atts)]
	# build function as sum for functions for each timepoint
	func_all <- function(h){ 0 }
	for(i in 1:ntimes){
		func_here <- function(h){ ((h*(1 - 2^(-D_prime[i]/(h*txnrate)))) - R_prime[i])^2 }
		func_all <- function(h){ func_all + func_here }
	}
	# initialize values for running optimizatin
	hold.row <- c(NA, NA)
	starth = 0
	if(sum(is.na(R_prime)) == ntimes){ return(hold.row) }
	# try running function
	try(fit.hold <- optim(starth, func_all), silent=T)
	try(hold.row <- c(fit.hold$par, fit.hold$value), silent=T)
  	return(hold.row)
}


###### running script #####

# get file with junction counts, Dprime, and Rprime
juncratio.data <- read.table(JUNCTIONFILE, sep="\t", header=T)

# separate out D' and R' values
D_primes <- juncratio.data[,grep('_Dprime', colnames(juncratio.data))]
R_primes <- juncratio.data[,grep('_Rprime', colnames(juncratio.data))]
combo_prime <- cbind(D_primes, R_primes)

# run model
sumsqfit.data <- t(apply(combo_prime, 1, sumsqequationsolve, TXNRATE))
# fitvalue = half live, yvalue = error around fit from optimization function
juncratio.data$fitvalue <- sumsqfit.data[,1]
names(juncratio.data)[names(juncratio.data) == 'fitvalue'] <- 'half.life'
juncratio.data$yvalue <- sumsqfit.data[,2]
names(juncratio.data)[names(juncratio.data) == 'yvalue'] <- 'fit.error'

##### Parse modeled introns to a set with optimal power to confidently detect half-lives
juncratio.data$result <- "success"
ratioinds <- grep('_ratio', colnames(juncratio.data))

# Remove those without intron-exon junction coverage in the earliest timepoint
## determined by ratio_T1 > 0 (if ie = 0, then ratio = 0)
#juncratio.data.parsed <- juncratio.data[which(juncratio.data[,ratioinds[1]] > 0),]
juncratio.data$result[which(juncratio.data[,ratioinds[1]] > 0)] <- "noIE_veryFast"

# Remove those without exon-exon junction reads in the last timepoint
## determined by ratio_Tn = inf
#juncratio.data.parsed <- juncratio.data[which(!is.infinite(juncratio.data[,ratioinds[length(ratioinds)]])),]
juncratio.data$result[which(is.infinite(juncratio.data[,ratioinds[length(ratioinds)]]))] <- "noEE_verySlow"

# Tag those where fit failed
#juncratio.data.parsed <- subset(juncratio.data, !is.na(fitvalue))
juncratio.data$result[which(is.na(juncratio.data$half.life))] <- "no_fit"

# Keep only those with coverage in at least n timepoints (where we used n=2)
## to be added?

##### write out final file with splicing half lives #####
#write.table(juncratio.data.parsed, file=JUNCTIONFILE, quote=F, sep="\t", row.names=F, col.names=T)
write.table(juncratio.data, file=JUNCTIONFILE, quote=F, sep="\t", row.names=F, col.names=T)
