MobilityFeatures <-
function(locations_df,
                            ACCURACY_LIM=51, ### meters GPS accuracy
                            ITRVL=10, ### seconds (data concatenation)
                            nreps=1, ### simulate missing data numer of times
                            tz="", ### time zone of data, defaults to current time zone
                            CENTERRAD=200, ### meters radius from significant locations considered
                            wtype="GLR",
                            spread_pars=c(10,1),
                            minpausedur=300,
                            minpausedist=60,
                            rad_fp=NULL,
                            wid_fp=NULL
){
  mobmatmiss=GPS2MobMat(locations_df,itrvl=ITRVL,accuracylim=ACCURACY_LIM,r=rad_fp,w=wid_fp)
  mobmat = GuessPause(mobmatmiss,mindur=minpausedur,r=minpausedist)
  obj=InitializeParams(mobmat)
  qOKmsg=MobmatQualityOK(mobmat,obj)
  if(qOKmsg!=""){
    cat(qOKmsg,"\n")
    return(NULL)
  }
  lsmf = list()
  lssigloc = list()
  for(repnum in 1:nreps){
    if(repnum==1){
      cat("Sim #: 1")
    }else if(repnum<=nreps-1){
      cat(paste(" ",repnum,sep=""))
    }else{
      cat(paste(" ",nreps,"\n",sep=""))
    }
    out3=SimulateMobilityGaps(mobmat,obj,wtype,spread_pars)
    IDundef=which(out3[,1]==3)
    if(length(IDundef)>0){
      out3=out3[-IDundef,]      
    }
    obj3=InitializeParams(out3)
    out_GMFM=GetMobilityFeaturesMat(out3,obj3,mobmatmiss,tz,CENTERRAD,ITRVL)
    lsmf[[repnum]]=out_GMFM[[1]]
    lssigloc[[repnum]]=out_GMFM[[2]]
  }
  cat("\n\n")
  if(length(lsmf)!=0){
    featavg = lsmf[[1]]
    if(nreps>1){
      for(i in 2:nreps){
        featavg=featavg+lsmf[[i]]
      }    
      featavg=featavg/nreps
    }    
  }else{
    featavg=NULL
  }
  return(list('mobmat'=mobmat,'mobmatmiss'=mobmatmiss,'featsims'=lsmf,'siglocsims'=lssigloc,'featavg'=featavg))
}
