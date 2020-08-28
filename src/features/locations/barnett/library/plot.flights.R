plot.flights <-
function(mat,xrang=NULL,yrang=NULL,diminch=6,add2plot=FALSE,addlegend=TRUE,outfile=NULL,title=NULL){
  #col24hour_v=c("#253494","#2c7fb8","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026","#7a0177","#ae017e","#dd3497","#f768a1","#fa9fb5","#fcc5c0","#edf8fb","#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c")
  col24hour_v=c(c("#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d"),rev(c("#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d")))  
  if(nrow(mat)==0){
    return(NULL)
  }
  if(add2plot){outfile=NULL}
  write2file=FALSE
  if(!is.null(outfile)){
    write2file=TRUE
  }
  if(write2file){
    pdf(outfile,width=diminch*11/10,height=diminch)
  }
  if(is.null(xrang)){
    xrang=plotlimits(mat)$xrang
  }
  if(is.null(yrang)){
    yrang=plotlimits(mat)$yrang
  }
  if(xrang[2]-xrang[1]>yrang[2]-yrang[1]){
    dif=((xrang[2]-xrang[1])-(yrang[2]-yrang[1]))
    yrang[2] = yrang[2]+dif/2
    yrang[1] = yrang[1]-dif/2
  }else{
    dif=((yrang[2]-yrang[1])-(xrang[2]-xrang[1]))
    xrang[2] = xrang[2]+dif/2
    xrang[1] = xrang[1]-dif/2    
  }
  xrang[1]=xrang[1]-(xrang[2]-xrang[1])/10
  if(!add2plot){
    if(!is.null(title)){
      par(mai=c(0,0,.4,0))      
      plot(NA,xlim=xrang
           ,ylim=yrang
           ,xaxt="n"
           ,yaxt="n"
           ,xlab=""
           ,ylab=""
           ,bty="n"
           ,main=title)  
    }else{
      par(mai=c(0,0,0,0))      
      plot(NA,xlim=xrang
           ,ylim=yrang
           ,xaxt="n"
           ,yaxt="n"
           ,xlab=""
           ,ylab=""
           ,bty="n"
           ,main="")        
    }
    xleg1=xrang[1]      
    xleg2=xleg1+(xrang[2]-xrang[1])/30
    yincr=(yrang[2]-yrang[1])/40
    legtext=c(" 6AM"," 12PM"," 6PM"," 12AM")
    if(addlegend){
      for(i in 1:24){
        polygon(c(xleg1,xleg1,xleg2,xleg2),c(yrang[1]+(i-1)*yincr,yrang[1]+i*yincr,yrang[1]+i*yincr,yrang[1]+(i-1)*yincr),col=col24hour_v[i])
        if(i%%6==0){
          text(xleg2,yrang[1]+i*yincr,legtext[floor(i/6)],adj=0,cex=.5)        
        }
      }
      points(xleg1,yrang[1]+26*yincr,cex=2,pch=16)
      text(xleg2,yrang[1]+26*yincr,">4 hrs",adj=0,cex=.5)
      points(xleg1,yrang[1]+28*yincr,cex=1,pch=16)
      text(xleg2,yrang[1]+28*yincr,"1 hr",adj=0,cex=.5)
      points(xleg1,yrang[1]+30*yincr,cex=.5,pch=16)
      text(xleg2,yrang[1]+30*yincr,"<30 mins",adj=0,cex=.5)
      #text(xleg1,yrang[1]+32.3*yincr,"Pause\nDuration",adj=0,cex=.5)
      legdist=10^floor(log10((yrang[2]-yrang[1])/5))
      lines(c(xrang[1],xrang[1]),c(yrang[2],yrang[2]-legdist))
      lines(c(xrang[1]-(xrang[2]-xrang[1])/120,xrang[1]+(xrang[2]-xrang[1])/120),c(yrang[2],yrang[2]))
      lines(c(xrang[1]-(xrang[2]-xrang[1])/120,xrang[1]+(xrang[2]-xrang[1])/120),c(yrang[2]-legdist,yrang[2]-legdist))
      if(log10(legdist)<3){
        text(xrang[1]+(xrang[2]-xrang[1])/50,yrang[2]-legdist/2,paste(as.character(legdist),"m"),cex=.5,srt=270,adj=c(.5,0))      
      }else{
        text(xrang[1]+(xrang[2]-xrang[1])/50,yrang[2]-legdist/2,paste(as.character(legdist/1000),"km"),cex=.5,srt=270,adj=c(.5,0))      
      }      
    }
  }
  for(i in 1:nrow(mat)){
    hour=as.numeric(strsplit(strsplit(as.character(as.POSIXct(mean(c(mat[i,4],mat[i,7]),na.rm=T),origin='1970-01-01'))," ")[[1]][2],":")[[1]][1])
    if(mat[i,1]==1){
      lines(c(mat[i,2],mat[i,5]),c(mat[i,3],mat[i,6]),col=col24hour_v[hour+1])
    }
    if(mat[i,1]==2){
      pwidth=max(.5,min(sqrt(2*(mat[i,7]-mat[i,4])/7200),2))
      points(mat[i,2],mat[i,3],pch=19,col=paste(col24hour_v[hour+1],"CC",sep=""),cex=pwidth)
    }
  }  
  if(write2file){
    dev.off()
  }
}
