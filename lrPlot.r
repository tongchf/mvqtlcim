Args <- commandArgs()
arglth <- length(Args)
if(arglth<7 || arglth>17 || arglth%%2==0){
cat("\nProgram: lrPlot.r (Plotting LR profiles in QTL mapping)\n")
cat("Version: 2016-05-10\n")
cat("Contact: Chunfa Tong <tongchf@njfu.edu.cn>\n\n")
cat("Usage: Rscript lrPlot.r -i qtlrstfile [options]\n\n")
cat("Options:\n")
cat("\t-p\tstr\tresult file for permutations\n")
cat("\t-w\tint\tjust taking the value of 0, 1 or 2 [0]\n") 
cat("\t\t\t0: manually giving the threshold value\n")
cat("\t\t\t1: genome-wide threshold determined by permutations\n")
cat("\t\t\t2: linkage-group-wide threshold determined by permutations\n")
cat("\t-s\tfloat\tsignificant level in(0.0001,0.2) [0.05]\n")
cat("\t-t\tstr\tthe plot format of 'pdf','png','jpg','tif' and 'bmp' [pdf]\n")
cat("\t-v\tfloat\tgiving the threshold value manually\n\n")
stop("Note: Parameter error!")
}

permurstfile=NULL
wglevel=0
siglevel=0.05
type="pdf"
threshold=NULL

for(i in seq(6,arglth,2)){
	if(Args[i]=="-i"){
		qtlrstfile=Args[i+1]
	}else if(Args[i]=="-p"){
		permurstfile=Args[i+1]
	}else if(Args[i]=="-w"){
		wglevel=as.numeric(Args[i+1])
	}else if(Args[i]=="-s"){
		siglevel=as.numeric(Args[i+1])
	}else if(Args[i]=="-t"){
		type=Args[i+1]
	}else if(Args[i]=="-v"){
		threshold=as.numeric(Args[i+1])
	}
}

#lrPlot <- function(qtlrstfile,permurstfile=NULL,
#			wglevel=0,siglevel=0.05,type="pdf",threshold=NULL){ 
# wglevel=0: manually setting the threshold when threshold must have a value
# wglevel=1: using whole genome threshold
# wglevel=2: using linkage group threshold
# type could be "pdf","bmp","png","jpg" or "tif"
  
  if(wglevel!=0 && wglevel!=1 && wglevel!=2){
	stop("The parameter wglevel value error!")
  }
  if(type!="pdf" && type!="bmp" && type!="png" && type!="jpg" && type!="tif"){
	stop("The parameter type value error!")
  }
  if(!file.exists(qtlrstfile)){
	stop("QTL mapping result file does not exist!")
  }
  dat <- read.table(qtlrstfile)

  if(is.null(permurstfile)){
	if(wglevel==1 || wglevel==2){
		stop("The parameter wglevel should be set to 0 when you did not provide the permutation result file!")
	}
  }else{
	if(wglevel==0){
		stop("The parameter wglevel should be set to 1 or 2 when you provide the permutation result file!")
	}
	if(!is.null(threshold)){
		stop("The parameter threshold should not be set when you provide the permutation result file!")
	}
  }

  if(siglevel<0.0001 && siglevel>0.2){
	stop("The parameter siglevel should be set in (0.0001,0.2)!")
  }

  ngs <- unique(dat[,1])
  ng <- length(ngs)
  gdst <- c()
  for(i in 1:ng){
  	gdst0 <- dat[dat$V1==ngs[i],4]
   	gdst <- c(gdst,min(gdst0)-0.01)
  }
  gdst <- c(gdst,max(gdst0)+0.01)

  map <- dat[,4]
  lrs <- dat[,7]
  vy <- ceiling(max(lrs))
  vx <- ceiling(max(map))
  if(vx < 1000){
	wdth <- 8
  }else{
  	wdth <- 8 + 2.5*(vx-1000)/1000
  }

 if(!is.null(permurstfile)){
  	permrst <- read.table(permurstfile)
  	npm <- ncol(permrst)-1
	if(npm<200){
		stop(paste("Note: The number of permutation times is only ",npm,", possibly too low!",sep=""))
	}
  	ulg <- unique(permrst[,1])
  	nlg <- length(ulg)

  	if(nlg!=ng || length(ngs)!=length(ulg)){
		stop("The mapping result file does not match the permutation result file!")
  	}
  	for(i in 1:ng){
		if(ngs[i]!=ulg[i]){
			stop("The mapping result file does not match the permutation result file!")
		}
  	}

  	lr <- c()
  	for(i in 1:npm){
		lr <- c(lr,max(permrst[,i+1]))
  	}
  	wgthd <- quantile(lr,1-siglevel)
  	lgthds <- c()
  	for(i in 1:nlg){
		lrsdat <- permrst[permrst$V1==ulg[i],]
		lrs0 <- c()
		for(j in 1:npm){
			lrs0 <- c(lrs0,max(lrsdat[,j+1]))
		}
		lgthds <- c(lgthds,quantile(lrs0,1-siglevel))
  	}
 }

  preoutfile = unlist(strsplit(qtlrstfile,'.',fixed=TRUE))[1]
  outfile = paste(preoutfile,"LRS","_W",wglevel,"_",siglevel,sep="")
  
  if(type=="pdf"){		
    	pdf(file=paste(outfile,".",type,sep=""),width=wdth,height=4)	   ## ploting pdf file
  }else if(type=="png"){
	png(file=paste(outfile,".",type,sep=""),width=wdth,height=4,units="in",res=600) ## ploting png file
  }else if(type=="jpg"){
	jpeg(file=paste(outfile,".",type,sep=""),width=wdth,height=4,units="in",res=600) ## ploting jpeg file
  }else if(type=="bmp"){
	bmp(file=paste(outfile,".",type,sep=""),width=wdth,height=4,units="in",res=300) ## ploting bmp file
  }else if(type=="tif"){
	tiff(file=paste(outfile,".",type,sep=""),width=wdth,height=4,units="in",res=600,
		compression="lzw+p") ## ploting tif file
  }else{
	print("Not support the plot type!")
  }
  
  plot(map,lrs,ylim=c(0,vy),type="l",lty=1,col="blue",xaxt ="n",
       yaxt="n",ylab="LR",xlab="Map position",las=1)
  itvx <- round(vx/600)*100
  itvy <- round(vy/25)*5
  axis(1,c(seq(0,vx,itvx)),c(seq(0,vx,itvx)));
  axis(2,c(seq(0,vy,itvy)),c(seq(0,vy,itvy)));
  
  for(i in 2:ng){
    abline(v=gdst[i],lty=3,col="gray")
    text((gdst[i-1]+gdst[i])/2,vy-0.5*itvy,i-1)
  }
  text((gdst[ng]+gdst[ng+1])/2,vy-0.5*itvy,ng)

  if(wglevel==1){
	abline(h=wgthd,lty=2)
  }else if(wglevel==2){
	segments(-1000,lgthds[1],x1=gdst[2],lty=2)
	for(i in 2:(ng-1)){
		segments(gdst[i],lgthds[i],x1=gdst[i+1],lty=2)	
	}
	segments(gdst[ng],lgthds[ng],x1=gdst[ng+1]+500,lty=2)
  }else{
	if(wglevel==0 && is.numeric(threshold)){
    		abline(h=threshold,lty=2)
	}
  }

  dev.off()   ## close the device

  # Summarizing QTL mapping results
  ugm <- unique(dat[,1:2])
  nugm <- nrow(ugm)
  outfile = paste(preoutfile,"_SigQtl_W",wglevel,"_",siglevel,".txt",sep="")
     
  if( wglevel==0 && is.numeric(threshold) ){
	head0 <- paste("# QTL mapping result for model ",dat[1,8], 
		" at the significant level of ",siglevel,"\n# with the whole-genome-wide threshold of ",
		threshold,"\n",sep="")
     write(head0,outfile)
	for(i in 1:nugm){
		dat0 <- dat[dat$V1==ugm[i,1] & dat$V2==ugm[i,2],]
		max0 <- max(dat0[,7])
		if(max0 >= threshold){
			rst0 <- dat0[dat0$V7==max0,]
     		write.table(rst0,outfile,append=TRUE, quote = FALSE, row.names =FALSE, col.names=FALSE)
		}	
	}
  }

 if( wglevel==1 ){
	head0 <- paste("# QTL mapping result for model ",dat[1,8], 
		" at the significant level of ",siglevel,"\n# with the whole-genome-wide threshold of ",
		wgthd,"\n",sep="")
     write(head0,outfile)
	for(i in 1:nugm){
		dat0 <- dat[dat$V1==ugm[i,1] & dat$V2==ugm[i,2],]
		max0 <- max(dat0[,7])
		if(max0 >= wgthd){
			rst0 <- dat0[dat0$V7==max0,]
     		write.table(rst0,outfile,append=TRUE, quote = FALSE, row.names =FALSE, col.names=FALSE)
		}	
	}
  }

  if( wglevel==2 ){
	head0 <- paste("# QTL mapping result for model ",dat[1,8], " at the significant level of ",
		siglevel,"\n# with the following linkage-group-wide thresholds:",sep="")
     write(head0,outfile)
	for(i in 1:ng){
		head0 <-paste("# LG",ngs[i],"\t",lgthds[i],sep="")
		write(head0,outfile,append=TRUE)
	}
	write(" ",outfile,append=TRUE)
	for(i in 1:nugm){
		dat0 <- dat[dat$V1==ugm[i,1] & dat$V2==ugm[i,2],]
		max0 <- max(dat0[,7])
		if(max0 >= lgthds[ugm[i,1]]){
			rst0 <- dat0[dat0$V7==max0,]
     		write.table(rst0,outfile,append=TRUE, quote = FALSE, row.names =FALSE, col.names=FALSE)
		}	
	}
  }
#}
