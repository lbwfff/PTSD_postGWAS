check_and_install <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

required_packages <- c("TeachingDemos", "MVMR", "tidyr",'tibble','dplyr',
                       'TwoSampleMR','coloc','ggplot2','ieugwasr','patchwork',
                       'geni.plots','plinkbinr') 

check_and_install(required_packages)

library("TeachingDemos")
library('MVMR')
library('tidyr')
library('tibble')
library('dplyr')
library('TwoSampleMR')
library('coloc')
library('ggplot2')
library('ieugwasr')
library('patchwork')
library('geni.plots')
library('plinkbinr')

# parameters for plot
genemove = 0.01; txt=1.1;  cex =1.3; lab=1.1; axis=1; top_cex=1.2;

GeneRowNum = function(GENELIST) {
  BP_THRESH = 0.03; MAX_ROW = 5
  # get the start and end position
  GENELIST = GENELIST[!duplicated(GENELIST$GENE),]
  START1 = as.numeric(GENELIST$GENESTART); END1 = as.numeric(GENELIST$GENEEND)
  STRLENGTH = nchar(as.character(GENELIST$GENE))
  MIDPOINT = (START1 + END1)/2
  START2 = MIDPOINT-STRLENGTH/250; END2 = MIDPOINT+STRLENGTH/250
  START = cbind(START1, START2); END = cbind(END1, END2);
  START = apply(START, 1, min); END = apply(END, 1, max)
  GENELIST = data.frame(GENELIST, START, END)
  GENELIST = GENELIST[order(as.numeric(GENELIST$END)),]
  START = as.numeric(GENELIST$START); END = as.numeric(GENELIST$END)
  # get the row index for each gene
  NBUF = dim(GENELIST)[1]
  ROWINDX = rep(1, NBUF)
  ROWEND = as.numeric(rep(0, MAX_ROW))
  MOVEFLAG = as.numeric(rep(0, NBUF))
  if(NBUF>1) {
    for( k in 2 : NBUF ) {
      ITERFLAG=FALSE
      if(START[k] < END[k-1]) {
        INDXBUF=ROWINDX[k-1]+1
      } else INDXBUF = 1
      if(INDXBUF>MAX_ROW) INDXBUF=1;
      REPTIME=0
      repeat{
        if( ROWEND[INDXBUF] > START[k] ) {
          ITERFLAG=FALSE
          INDXBUF=INDXBUF+1
          if(INDXBUF>MAX_ROW) INDXBUF = 1
        } else {
          ITERFLAG=TRUE
        }
        if(ITERFLAG) break;
        REPTIME = REPTIME+1
        if(REPTIME==MAX_ROW) break;
      }
      ROWINDX[k]=INDXBUF;
      
      if( (abs(ROWEND[ROWINDX[k]]-START[k]) < BP_THRESH)
          | ((ROWEND[ROWINDX[k]]-START[k])>0) ) {
        MOVEFLAG[k] = 1
        SNBUF = tail(which(ROWINDX[c(1:k)]==ROWINDX[k]), n=2)[1]
        MOVEFLAG[SNBUF] = MOVEFLAG[SNBUF] - 1
      }
      if(ROWEND[ROWINDX[k]]<END[k]) {
        ROWEND[ROWINDX[k]] = END[k]  }
    }
  }
  GENEROW = data.frame(as.character(GENELIST$GENE),
                       as.character(GENELIST$ORIENTATION),
                       as.numeric(GENELIST$GENESTART),
                       as.numeric(GENELIST$GENEEND),
                       ROWINDX, MOVEFLAG)
  colnames(GENEROW) = c("GENE", "ORIENTATION", "START", "END", "ROW", "MOVEFLAG")
  return(GENEROW)
}

plot_probe = function(probeinfobuf, k, colplot, x.min, x.max, y.min, y.max,pchbuf,heidi) {
  xcenter = as.numeric(probeinfobuf[k,3])
  pvalbuf = as.numeric(probeinfobuf[k,8])
  strbuf = probeinfobuf[k,1]
  par(new=TRUE)
  if(heidi==TRUE) {
    plot(xcenter, pvalbuf, ylim=c(y.min,y.max),  xlim=c(x.min,x.max),cex.axis=axis,
         xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
  } else {
    plot(xcenter, pvalbuf, ylim=c(y.min,y.max),  xlim=c(x.min,x.max),cex.axis=axis,
         xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
  }
}

ReadSMRData = function(plotfile)
{
  SMRData = list();
  key=c("$probe","$SNP","$GWAS","$eQTL");
  skiplines=0;
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(keywords[1]!=key[1])
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  nprobes=as.numeric(keywords[2]);
  SMRData$probeID=keywords[3];
  
  
  skiplines=skiplines+1;
  SMRData$SMR=read.table(plotfile, header=F, nrows=nprobes, skip=skiplines);
  skiplines=skiplines+nprobes;
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(keywords[1]!=key[2])
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  nrs=as.numeric(keywords[2]);
  skiplines=skiplines+1;
  SMRData$SNP=read.table(plotfile, header=F, nrows=nrs, skip=skiplines);
  skiplines=skiplines+nrs;
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(keywords[1]!=key[3])
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  ngwas=as.numeric(keywords[2]);
  skiplines=skiplines+1;
  SMRData$GWAS=read.table(plotfile, header=F, nrows=ngwas, skip=skiplines);
  skiplines=skiplines+ngwas;
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(keywords[1]!=key[4])
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  neqtl=as.numeric(keywords[2]);
  skiplines=skiplines+1;
  
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  prbname=keywords[1];
  neqtlsnp=as.numeric(keywords[2]);
  skiplines=skiplines+1;
  SMRData$eQTL=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
  SMRData$eQTL=cbind(prbname,SMRData$eQTL)
  skiplines=skiplines+neqtlsnp;
  if(neqtl>1)
  {
    for(i in 2:neqtl)
    {
      keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
      prbname=keywords[1];
      neqtlsnp=as.numeric(keywords[2]);
      skiplines=skiplines+1;
      raweQTLtmp=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
      raweQTLtmp=cbind(prbname,raweQTLtmp);
      SMRData$eQTL=rbind(SMRData$eQTL,raweQTLtmp);
      skiplines=skiplines+neqtlsnp;
    }
  }
  
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(length(keywords)>0)
  {
    if(keywords[1]!="$Gene")
    {
      print("ERROR: plot file is not correct!");
      quit();
    }
    ngenes=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$Gene=read.table(plotfile, header=F, nrows=ngenes, skip=skiplines);
  }
  return(SMRData)
}

SMRLocusPlot = function(data=SMRData, probeNEARBY=NULL,smr_thresh=NULL, smr_thresh_plot=NULL, heidi_thresh=NULL, plotWindow=NULL,pointsize=20,max_anno_probe=16,anno_selfdef=TRUE)
{
  
  cex_coeff=3/4 * pointsize/15;
  if(length(smr_thresh)==0){
    print("ERROR: please specify the threshold of SMR test!");
    quit();
  }
  if(length(heidi_thresh)==0){
    print("ERROR: please specify the threshold of HEIDI test!");
    quit();
  }
  if(length(plotWindow)==0){
    print("ERROR: please specify the plot window size!");
    quit();
  }
  if(length(which(is.na(data$SMR[,3])))>0)
  {
    print("ERROR: Some probes' physical positon is missing!");
    quit();
  }
  idx=match(data$probeID,data$SMR[,1]);
  if(length(idx)==0){
    print("ERROR: Plot file is not generated correctly, can't find target probe!");
    quit();
  }
  if(length(smr_thresh_plot)==0){
    smr_thresh_plot=smr_thresh;
  }
  cis_start=data$SMR[idx,3]-plotWindow*1000;
  if(cis_start<0) cis_start=0
  cis_end=data$SMR[idx,3]+plotWindow*1000;
  idx=which(data$SMR[,3]>=cis_start & data$SMR[,3]<=cis_end)
  data$SMR=data$SMR[idx,]
  idx=match(data$GWAS[,1],data$SNP[,1])
  tmpsnpbp=data$SNP[idx,3]
  idx=which(tmpsnpbp>=cis_start &tmpsnpbp<=cis_end)
  data$GWAS=data$GWAS[idx,]
  idx=match(data$eQTL[,2],data$SNP[,1])
  tmpsnpbp=data$SNP[idx,3]
  idx=which(tmpsnpbp>=cis_start &tmpsnpbp<=cis_end)
  data$eQTL=data$eQTL[idx,]
  
  if(!is.null(data$Gene))
  {
    idx=which(data$Gene[,2]>=cis_start & data$Gene[,3]<=cis_end )
    data$Gene=data$Gene[idx,]
  }
  
  #start to plot
  smrindx = which(data$SMR[,8] <= smr_thresh_plot)
  #heidiindx = which((data$SMR[,8] <= smr_thresh_plot) & (data$SMR[,9] >= heidi_thresh_plot))
  smrprobes = NULL; heidiprobes = NULL;
  if(length(smrindx)>0) { smrprobes =  as.character(data$SMR[smrindx,1]) }
  #if(length(heidiindx)>0) { heidiprobes = as.character(data$SMR[heidiindx,1]) }
  
  smrindx_bonferr = which(data$SMR[,8] <= smr_thresh)
  heidiindx_strengent = which((data$SMR[,9] >= heidi_thresh))
  smrprobes_red = NA; heidiprobes_solid = NA;
  if(length(smrindx_bonferr)>0) { smrprobes_red =  as.character(data$SMR[smrindx_bonferr,1]) }
  if(length(heidiindx_strengent)>0) { heidiprobes_solid = as.character(data$SMR[heidiindx_strengent,1]) }
  
  if(length(probeNEARBY)>0)
  {
    idx=match(probeNEARBY,data$SMR[,1])
    idxx=which(is.na(idx))
    if(length(idxx)>0)
    {
      for(ii in 1:length(idxx)) {
        print(paste("WARNING: cann't find probe ",probeNEARBY[idxx[ii]], " in plot region.",sep=""))
      }
      probeNEARBY=probeNEARBY[-idxx]
    }
    
  }
  probePLOT=smrprobes #draw the eQTL of all the probes that passed smr_thresh_plot
  probePLOT=unique(c(data$probeID,probePLOT,probeNEARBY)) # draw the target probe anyway
  nprobePLOT = length(probePLOT)
  
  idx=which(is.na(data$GWAS[,2]) | is.na(data$GWAS[,3]))
  if(length(idx)>0) data$GWAS=data$GWAS[-idx,]
  pZY=-log10(pchisq((data$GWAS[,2]/data$GWAS[,3])^2,1,lower.tail=F))
  
  idx=match(data$probeID,data$SMR[,1]);
  if(length(idx)>0){
    chrPLOT = data$SMR[idx,2]
  }else{
    print("ERROR: Plot file is not generated correctly, please report this bug!");
    quit();
  }
  idx=which(is.na(data$SMR[,8]) )
  if(length(idx)>0) {
    probeINFO=data$SMR[-idx,];
  }else{
    probeINFO=data$SMR;
  }
  idx=which(is.na(probeINFO[,5]) | is.na(probeINFO[,6]));
  idx2=which(is.na(probeINFO[,3]));
  if(length(intersect(idx,idx2))>0)
  {
    print("ERROR: Some probes' physical positon is missing!");
    quit();
  }
  probeINFO[idx,5]=probeINFO[idx,3]-7500;
  probeINFO[idx,6]=probeINFO[idx,3]+7500;
  probeINFO[,8]=-log10(probeINFO[,8]);
  probeINFO[,3]=probeINFO[,3]/1e6;
  probeINFO[,5]=probeINFO[,5]/1e6;
  probeINFO[,6]=probeINFO[,6]/1e6;
  pXY=probeINFO[,8];
  yMAX = ceiling(max(c(pZY, pXY), na.rm=T)) + 1;
  if(is.null(data$Gene))
  {
    glist=cbind(probeINFO[,2],probeINFO[,5:6],as.character(probeINFO[,4]),probeINFO[,7]);
  } else {
    glist=data$Gene;
    glist[,2]=glist[,2]/1e6;
    glist[,3]=glist[,3]/1e6;
  }
  colnames(glist)=c("CHR", "GENESTART",  "GENEEND",   "GENE", "ORIENTATION");
  idx=which(is.na(glist[,2]) | is.na(glist[,3]));
  if(length(idx>0)) glist=glist[-idx,];
  generow = GeneRowNum(glist);
  num_row = max(as.numeric(generow$ROW));
  offset_map = ceiling(yMAX);
  offset_probe = yMAX / 2.5;
  num_probe = nprobePLOT
  offset_eqtl = ceiling(yMAX / 2.5) + 0.5;
  dev_axis = 0.1*yMAX;
  if(dev_axis<1.5) dev_axis = 1.5;
  yaxis.min = -offset_map - offset_eqtl*num_probe - dev_axis*(num_probe+1);
  yaxis.max = yMAX + ceiling(offset_probe) + 1;
  # scales of x-axis
  idx=match(data$GWAS[,1],data$SNP[,1]);
  gwasBP = as.numeric(data$SNP[idx,3])/1e6;
  #min.pos = min(gwasBP);
  #max.pos = max(gwasBP);
  min.pos = cis_start/1e6
  max.pos = cis_end/1e6
  start = min(as.numeric(glist[,2]));
  end = max(as.numeric(glist[,3]));
  bp = c(min.pos, max.pos, start, end);
  xmin = min(bp, na.rm=T) - 0.001;  xmax = max(bp, na.rm=T) +0.001;
  xmax=xmax+(xmax-xmin)*0.1 #extend
  ylab = expression(paste("-", log[10], "(", italic(P), " GWAS or SMR)", sep=""));
  xlab = paste("Chromosome", chrPLOT, "Mb");
  # plot GWAS p value
  par(mar=c(5,5,3,2), xpd=TRUE)
  plot(gwasBP, pZY, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max),
       ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
       xlim=c(xmin, xmax), pch=20, col="gray68");
  
  # y1 axis
  devbuf1 = yMAX/4
  axis(2, at=seq(0,yMAX,devbuf1), labels=round(seq(0,yMAX,devbuf1),0), las=1, cex.axis=axis);
  mtext(ylab, side=2, line=3, at=(yMAX*2/3), cex=cex_coeff);
  eqtl.lab = expression(paste("-", log[10], "(", italic(P), " eQTL)", sep=""));
  axis.start = 0; axis.down = offset_eqtl + dev_axis;
  for( k in 1 : nprobePLOT ) {
    axis.start = axis.start - axis.down
    eqtlinfobuf = data$eQTL[which(data$eQTL[,1]==probePLOT[k]),]
    if(dim(eqtlinfobuf)[1]==0) next;
    pvalbuf=-log10(pchisq((eqtlinfobuf[,3]/eqtlinfobuf[,4])^2,1,lower.tail=F));
    pvalbuf[which(is.infinite(pvalbuf))]=1e-300;
    if(length(which(smrprobes_red==probePLOT[k]))==0) {
      col_eqtl = "navy"
    } else col_eqtl = "maroon"
    eqtl.min = 0; eqtl.max = ceiling(max(pvalbuf))
    eqtl.max =ceiling(eqtl.max *1.25) #extend
    pvalbuf = pvalbuf/eqtl.max * offset_eqtl + axis.start
    idx=match(eqtlinfobuf[,2],data$SNP[,1]);
    eqtlbp = as.numeric(data$SNP[idx,3])/1e6;
    probegene = unique(as.character(data$SMR[which(data$SMR[,1]==probePLOT[k]),4]))
    par(new=TRUE)
    pchbuf = 4;
    #if(k%%2==0) pchbuf = 20;
    plot(eqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
         ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
    # annotate the eQTLs
    text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=probePLOT[k], geneid=probegene)),col="black", cex=1, adj=0)
    # axis
    devbuf1 = offset_eqtl/3; devbuf2 = eqtl.max/3
    axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
         labels=round(seq(0,eqtl.max,devbuf2),0),
         las=1, cex.axis=axis)
    # add separator line
    segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
             col="dim grey", lty="24", lwd=1)
  }
  #ypos = (axis.start - dev_axis)/2
  ypos = (axis.start - dev_axis)*2/3
  mtext(eqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff)
  
  # plot p value of bTG
  # all the probes
  num_gene = dim(generow)[1]
  dist = offset_map/num_row
  for( k in 1 : num_row ) {
    generowbuf = generow[which(as.numeric(generow[,5])==k),]
    xstart = as.numeric(generowbuf[,3])
    xend = as.numeric(generowbuf[,4])
    snbuf = which(xend-xstart< 1e-3)
    if(length(snbuf)>0) {
      xstart[snbuf] = xstart[snbuf] - 0.0025
      xend[snbuf] = xend[snbuf] + 0.0025
    }
    xcenter = (xstart+xend)/2
    xcenter = spread.labs(xcenter, mindiff=0.01, maxiter=1000, min = xmin, max = xmax)
    num_genebuf = dim(generowbuf)[1]
    for( l in 1 : num_genebuf ) {
      ofs=0.3
      if(l%%2==0) ofs=-0.8
      m = num_row - k
      ypos = m*dist + yaxis.min
      code = 1
      if(generowbuf[l,2]=="+") code = 2;
      arrows(xstart[l], ypos, xend[l], ypos, code=code, length=0.07, ylim=c(yaxis.min,yaxis.max),
             col=colors()[75], lwd=1)
      movebuf = as.numeric(generowbuf[l,6])*genemove
      text(xcenter[l]+movebuf, ypos,label=substitute(italic(genename), list(genename=as.character(generowbuf[l,1]))), pos=3, offset=ofs, col="black", cex=0.9)
    }
  }
  
  # plot the probes
  probeINFO=probeINFO[order(probeINFO[,8],decreasing = TRUE),];
  nprobeINFO=dim(probeINFO)[1];
  if(nprobeINFO>max_anno_probe){
    probeINFO=probeINFO[c(1:max_anno_probe),]
    nprobeINFO=dim(probeINFO)[1];
  }
  if(anno_selfdef) probeINFO = probeINFO[order(probeINFO[[2]], probeINFO[[3]]), ]
  #probeINFO=probeINFO[order(probeINFO[2],probeINFO[3]),] ####20170217
  xcenter = as.numeric(probeINFO[,3])
  xcbuf = xcenter
  ####20170217####
  if(anno_selfdef)
  {
    reginlength=(xmax-(xmax-xmin)*0.15)-xmin
    leftspot=xmin+reginlength/20
    rightspot=(xmax-(xmax-xmin)*0.15)-reginlength/20
    itvl=(rightspot-leftspot)/dim(probeINFO)[1]
    if(dim(probeINFO)[1]==1) {
      xcenter=as.numeric(probeINFO[,3])
    } else {
      xcenter=leftspot+itvl/2
      for( k in 2:dim(probeINFO)[1]) xcenter=c(xcenter,leftspot+k*itvl)
    }
    
  } else {
    xcenter = spread.labs(xcenter[1:nprobeINFO], mindiff=0.08, maxiter=1000, min = xmin, max = xmax-1)
  }
  # adjust the line position
  
  adjflag = rep(0, nprobeINFO)
  if(nprobeINFO>1) {
    dbuf = c(0, xcbuf[1:(nprobeINFO-1)])
    mflag = as.numeric(abs(xcbuf[1:(nprobeINFO)] - dbuf) < 0.01)
    adjflag = as.numeric( mflag | c(mflag[2:nprobeINFO],0) )
  }
  
  for( k in 1 : nprobeINFO)  {
    hitflag=FALSE
    if(length(which(heidiprobes_solid==probeINFO[k,1]))>0 & length(which(smrprobes_red==probeINFO[k,1]))>0) {
      hitflag=TRUE
      colplot = "maroon"; colfont=2; pchbuf=23;
    } else if(length(which(smrprobes_red==probeINFO[k,1]))>0) {
      colplot = "maroon"; colfont=2; pchbuf=5
      #} else if (length(which(heidiprobes_solid==probeINFO[k,1]))>0) {
      #hitflag=TRUE
      # colplot = "navy"; colfont=1; pchbuf=23
    } else {
      colplot = "navy"; colfont=1; pchbuf=5
    }
    if( as.numeric(probeINFO[k,8]) < 0 ) {
      colplot = "black"; colfont=1;
    }
    # plot p value of bxy
    plot_probe(probeINFO, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
    # annotate the probes
    if(k<=max_anno_probe)
    {
      ypos = 1.02*yMAX
      strbuf =
        text(xcenter[k], ypos,
             labels=substitute(paste(probeid, " (", italic(genename), ")", sep=""),
                               list(probeid=as.character(probeINFO[k,1]),
                                    genename=as.character(probeINFO[k,4]))),
             ylim=c(yaxis.min, yaxis.max),
             srt=30, col=colplot, font=colfont, cex=1, adj=0)
      # plot the lines
      # 1st step
      xstart = xcbuf[k]
      ystart = as.numeric(probeINFO[k,8]); yend = yMAX*(1-1/20);
      if( nprobeINFO > 1 ) {
        if(adjflag[k]==1) {
          xstart = (xcbuf[k] + xcenter[k])/2
          segments(xcbuf[k], ystart, xstart, ystart, col=colplot, lwd=axis, lty=3)
        }
      }
      segments(xstart, ystart, xstart, yend, col=colplot, lwd=axis, lty=3)
      # 2nd step
      xend = xcenter[k]; ystart = yMAX*(1-1/20); yend = yMAX*1.01;
      segments(xstart, ystart, xend, yend, col=colplot, lwd=axis, lty=3)
    }
  }
  # plot the threshold
  # SMR threshold
  ybuf = -log10(as.numeric(smr_thresh)); dev_anno = yMAX/9;
  strbuf = paste("pSMR = ",smr_thresh, sep="")
  segments(xmin, ybuf, xmax, ybuf, col="maroon", lty=2, lwd=1);
  text(xmax, ybuf+dev_anno, labels=strbuf, adj=1, col="maroon", cex=axis,font=3);
}

tidy_mvmr_output <- function(mvmr_res) {
  #  tidy up MVMR returned output
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    dplyr::rename(b='Estimate',  
                  se="Std. Error",
                  pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=4) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 4),
           pval=as.numeric(pval))
}  

make_mvmr_input <- function(exposure_dat, outcome.id.mrbase=NULL, outcome.data=NULL){
  # provide exposure_dat created in the same way as for TwoSampleMR 
  # also specify the outcome argument [only ONE!] (MR-base ID or full gwas data in .outcome format)
  
  # extract SNPs for both exposures from outcome dataset
  # (for the selected option mr.base or local outcome data)
  if (!is.null(outcome.id.mrbase)) {
    # if mrbase.id is provided
    outcome_dat <- extract_outcome_data(snps = unique(exposure_dat$SNP),
                                        outcomes = outcome.id.mrbase)
  } else if (!is.null(outcome.data)){
    # if outcome df is provided
    outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  }
  
  # harmonize datasets
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  # Create variables for the analysis 
  
  ### works for many exposures
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  # add beta/se names
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}

MR_and_coloc <- function (SMRfile,trait_name,plotpath,
                          n_qtl,n_gwas,g_version,all_outcome,
                          bfile) {
  
  table<-list()
  
  ####################
  #SMR plot

  SMRData <- ReadSMRData(SMRfile) 
  
  probe<-SMRData[["probeID"]]
  
  path<-paste0(plotpath,'/',trait_name,'_',SMRData[["SMR"]]$V4[SMRData[["SMR"]]$V1==probe],'.SMRLocusPlot.pdf')

  pdf(path,width = 12,height = 8)
  SMRLocusPlot(data=SMRData, smr_thresh=0.01, heidi_thresh=0.05, plotWindow=600, max_anno_probe=16)
  dev.off()
  
  ##################
  #MR
  
  gwas<-SMRData[["GWAS"]]
  qtl<-SMRData[["eQTL"]]
  snp<-SMRData[["SNP"]]
  gwas<-merge(gwas,snp,by='V1')
  gwas$p<-(pchisq((gwas[,2]/gwas[,3])^2,1,lower.tail=F))
  
  qtl<-merge(qtl,snp,by='V1')
  qtl<-qtl[qtl$prbname==probe,]
  qtl$p<-(pchisq((qtl[,3]/qtl[,4])^2,1,lower.tail=F))
  
  exp_dat <- format_data(
    data.frame(qtl),type = "exposure",
    snp_col = "V1",beta_col = "V2.x",se_col = "V3.x",
    effect_allele_col = "V4.y",
    other_allele_col = "V5",
    pval_col = "p",gene_col = "gene_id",
    chr_col = "V2.y",pos_col = "V3.y")
  
  exp_dat$id.exposure<-SMRData[["SMR"]]$V4[SMRData[["SMR"]]$V1==probe]

  gwas<-merge(gwas,all_outcome,by.x=c('V1'),by.y=c('rsids'))
  
  outcome <- format_data(
    data.frame(gwas),type="outcome",
    snps=exp_dat$SNP,snp_col = "V1",
    beta_col = "V2.x",se_col = "V3.x",
    effect_allele_col = "V4",
    other_allele_col = "V5",
    pval_col = "p",eaf_col = 'af_alt')
  
  outcome$id.outcome<-trait_name
  
  coloc<-merge(exp_dat,outcome,by=c('SNP'))
  
  exp_dat<-exp_dat[exp_dat$pval.exposure<(10^-6),]
  
  clump <- ieugwasr::ld_clump(clump_r2=0.2,
                              dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id='test'),
                              plink_bin = plinkbinr::get_plink_exe(),
                              bfile = bfile)
  
  exp_dat<-exp_dat[exp_dat$SNP %in% clump$rsid,]
  
  MRinput <- harmonise_data(
    exposure_dat = exp_dat, 
    outcome_dat = outcome
  )
  MRinput$samplesize.outcome<-n_gwas
  
  df_MR_hetero <- mr_heterogeneity(MRinput)  
  df_MR_pleio <- mr_pleiotropy_test(MRinput) 
  
  table[['hetero']]<-df_MR_hetero
  table[['plei']]<-df_MR_pleio
  
  if (df_MR_hetero$Q_pval[nrow(df_MR_hetero)]<0.05) {
    res=mr(MRinput, method_list = c("mr_ivw_mre",
                                    "mr_egger_regression","mr_weighted_median", 
                                    'mr_simple_mode','mr_weighted_mode'))} else{
                                      res=mr(MRinput, method_list = c("mr_ivw_fe",
                                                                      "mr_egger_regression","mr_weighted_median", 
                                                                      'mr_simple_mode','mr_weighted_mode'
                                      ))}
  
  R2 <- MRinput$beta.exposure^2/(MRinput$beta.exposure^2 + (MRinput$se.exposure^2)*n_qtl)
  Fs <- R2*(n_qtl - nrow(MRinput) - 1)/nrow(MRinput)*(1 - R2)
  table[['F-stats']]<-Fs
  
  OR <-generate_odds_ratios(res)
  table[['OR']]<-OR
  
  plot<-list()
  
  res_loo <- mr_leaveoneout(MRinput)
  p<-mr_leaveoneout_plot(res_loo)[[1]]
  plot[[1]]<-p+scale_color_brewer(palette = 'Set2')+
    xlab('MR leave−one−out sensitivity analysis')+
    theme_classic(base_size = 16)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.background = element_blank(),
          legend.position = 'none',
          # legend.title = element_blank(),
          # legend.position = c(0.1,0.9),
          legend.text = element_text(size = 18))+
    theme(aspect.ratio=1.5)
  
  p<-mr_scatter_plot(res, MRinput)[[1]]
  p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
  p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
  p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
  plot[[2]]<-p+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 14)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.background = element_blank(),
          legend.position = 'top',
          # legend.title = element_blank(),
          # legend.position = c(0.1,0.9),
          legend.text = element_text(size = 12))+ 
    theme(aspect.ratio=1)+
    scale_x_continuous(limits = c(0, NA), expand = c(0, 0))
  
  res_single <- mr_singlesnp(MRinput)
  
  p<-mr_forest_plot(res_single)[[1]]
  p[["data"]]<-p[["data"]][-nrow(p[["data"]]),] 
  p[["layers"]][[4]]<-NULL
  
  plot[[3]]<-p+scale_color_brewer(palette = 'Set2')+
    geom_hline(yintercept='', linetype='dotted', col = 'black')+
    xlab('MR effect size')+
    theme_classic(base_size = 16)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.background = element_blank(),
          legend.position = 'none',
          # legend.title = element_blank(),
          # legend.position = c(0.1,0.9),
          legend.text = element_text(size = 18)) +
    theme(aspect.ratio=1.5)
  
  p<-mr_funnel_plot(res_single)[[1]]
  p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
  plot[[4]]<-p+scale_color_brewer(palette = 'Set1')+
    xlab('MR leave−one−out sensitivity analysis')+
    theme_classic(base_size = 16)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.background = element_blank(),
          legend.position = 'top',
          # legend.title = element_blank(),
          # legend.position = c(0.1,0.9),
          legend.text = element_text(size = 18))+
    theme(aspect.ratio=1)
  
  
  path<-paste0(plotpath,'/',trait_name,'_',SMRData[["SMR"]]$V4[SMRData[["SMR"]]$V1==probe],'.Twosample_MR_Plot.pdf')
  
  pdf(path,width = 12,height = 8)
  print(wrap_plots(plot,nrow=2) )
  dev.off()

  ####################
  #coloc
  
  coloc<-merge(coloc,all_outcome,by.x=c('SNP'),by.y=c('rsids'))
  coloc<-coloc[!duplicated(coloc$SNP),]
  
  D2<-list(type='quant',beta=coloc$beta.exposure,snp=coloc$SNP,
           varbeta=c(coloc$se.exposure)^2,sdY=1,N=n_qtl,
           position=coloc$pos.exposure)
  
  D1<-list(type='cc',beta=coloc$beta.outcome,varbeta=(coloc$se.outcome)^2,
           snp=coloc$SNP,  
           position=coloc$pos.exposure,N=n_gwas,MAF=coloc$af_alt)
  
  check_dataset(D2)
  check_dataset(D1)
  
  my.coloc <- coloc.abf(dataset1=D1,
                        dataset2=D2)
  
  highlights<-subset(my.coloc$results,SNP.PP.H4>0.6)$snp
  
  print(my.coloc$summary)
  
  table[['coloc']]<-my.coloc
  
  path<-paste0(plotpath,'/',trait_name,'_',SMRData[["SMR"]]$V4[SMRData[["SMR"]]$V1==probe],'.sensitivity.pdf')
  pdf(path,width = 12,height = 8)
  sensitivity(my.coloc,rule="H4 > 0.6") 
  dev.off()

  assoc<-coloc[,c('SNP','chr.exposure','pos.exposure','pval.outcome','pval.exposure')]
  colnames(assoc)<-c('marker','chr','pos','pvalue_1','pvalue_2')
  
  corr<-ieugwasr::ld_matrix(assoc$marker, with_alleles = F, 
                            plink_bin = plinkbinr::get_plink_exe(),
                            bfile = bfile)
  
  assoc<-assoc[match(rownames(corr),assoc$marker),]
  assoc$pos<-as.integer(assoc$pos)
  
  path<-paste0(plotpath,'/',trait_name,'_',SMRData[["SMR"]]$V4[SMRData[["SMR"]]$V1==probe],'.fig_region_stack.pdf')
  
  pdf(path,width = 12,height = 8)
  print(fig_region_stack(
    data = assoc,
    traits = c("PTSD", SMRData[["SMR"]]$V4[SMRData[["SMR"]]$V1==probe]),
    corr = corr,
    build = g_version,
    highlights = highlights,
    title_center = TRUE))
  dev.off()
  
  return(table)

}





