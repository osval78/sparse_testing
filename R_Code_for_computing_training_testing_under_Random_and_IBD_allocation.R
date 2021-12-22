rm(list=ls(all=TRUE))
library(crossdes)
library(dplyr)
#nP: No. of partitions
#nBa: No of Blocks to form in the IBD
#2 cases: Random, IBD
#To ensure that lambda>=1, then nBa*(k-1)/trt>2 ==> k>trt/nBa+2

PT_Tr_IBD_f<-function(dat,GID='GID',Block_Name='Loc',nP=10,nBa=2,ss=1e4)
{
  GIDs =  unique(dat[,GID])
  Block_Names =  unique(dat[,Block_Name])
  #by(dat$Gid,dat$Loc,function(x)length(unique(x)))
  set.seed(ss)
  PT_Tr =  data.frame()
  for(p in 1:nP)
  {
    print(paste("Particion ",p," de ",nP,sep=""))
    print("PT en proceso")
    k=ceiling(length(GIDs)*nBa/length(Block_Names))
    PT = find.BIB(trt=length(GIDs),b=length(Block_Names), k=k)
    print("PT terminado")
    PT2 = data.frame(Block=rep(1:length(Block_Names),each=dim(PT)[2]),
                     GID = c(t(PT)))
    print("PT2 terminado")
    PT2[,Block_Name] = Block_Names[PT2[,'Block']]
    PT2[,GID] = GIDs[PT2$GID]
    head(PT2)
    PT_Tr =  rbind(PT_Tr,data.frame(PT=p,Block_Name= PT2[,Block_Name],GID=PT2[,GID]))
  }
  No_BLocks_by_GID = PT_Tr%>%group_by(PT,GID)%>%select(Block_Name)%>%summarise(No_Blocks=length(unique(Block_Name)))
  No_BLocks_by_GID =  data.frame(No_BLocks_by_GID)
  PT_Tr_ls =  list(PT_Tr=PT_Tr, No_BLocks_by_GID=No_BLocks_by_GID)
  PT_Tr_ls
}

PT_tst_IBD_ls_f<-function(dat,GID='GID',Block_Name='Loc',nP=10,nBa=2,ss=1e4)
{
  print("working on training")
  PT_Tr_ls= PT_Tr_IBD_f(dat,GID,Block_Name,nP,nBa)
  #head(PT_Tr_ls$PT_Tr)
  #No_BLocks_by_GID = PT_Tr_ls$No_BLocks_by_GID
  #barplot(No_BLocks_by_GID$No_Blocks[No_BLocks_by_GID$PT==1])
  Pos_tr =  match(paste(PT_Tr_ls$PT_Tr$Block_Name,PT_Tr_ls$PT_Tr$GID,sep='_'),
                  paste(dat[,Block_Name],dat[,GID],sep='_'))
  PT_Tr_ls$PT_Tr$Pos_tr  =  Pos_tr
  PT_tr = PT_Tr_ls$PT_Tr
  head(PT_tr)
  PT_tst_ls =  vector('list',nP)
  print("working on testing")
  for(p in 1:nP)
  {
    print(paste("Particion ",p," de ",nP,sep=""))
    PT_tr_p =  PT_tr[PT_tr$PT==p,]
    PT_tst_ls[[p]] = setdiff(1:dim(dat)[1],PT_tr_p$Pos_tr)
  }
  names(PT_tst_ls) = paste0('PT',1:nP)
  PT_tst_ls
}

#Random assignation of each GID to each Block
RIBD_f0<-function(dat,GID='GID',Block_Name='Loc',nBa=2)
{
  dat[,GID] =  as.character(dat[,GID])
  GIDs =  unique(dat[,GID])
  Block_Names =  unique(dat[,Block_Name])
  nB =  length(Block_Names)
  k=ceiling(length(GIDs)*nBa/length(Block_Names))
  
  RIBD =  matrix(NA,nc=k,nr=nB)
  GIDs_R = GIDs
  GIDs_W = data.frame(GID=GIDs,W=nBa)
  for(j in 1:nB)
  {
    if(sum(GIDs_W$W>0)>=k)
    {
      RIBD[j,]  = sample(GIDs_W$GID,k,prob=GIDs_W$W)
      Tb_R = table(RIBD)
      Tb_R = data.frame(GID = names(Tb_R),Freq = c(Tb_R))
      if(dim(Tb_R)[1]>0)
      {
        GIDs_W$W[match(Tb_R$GID,GIDs_W$GID)] =nBa-Tb_R$Freq
      }
    } else
    {
      if(nB-j==0)
      {
        RIBD[j,]  = sample(GIDs_W$GID,k)
      }
      else
      {
        j=nB+1
        break
      }
      
    }
  }
  RIBD = data.frame(Block=Block_Names,RIBD)
  ifelse(j==nB+1,return("Not enough GIDs"), return(RIBD))
}

PT_Tr_RIBD_f<-function(dat,GID='GID',Block_Name='Loc',nP=10,nBa=2,ss=1e4)
{
  set.seed(ss)
  RIBD = data.frame()
  for(p in 1:nP)
  {
    RIBD_p = RIBD_f0(dat,GID,Block_Name,nBa)
    RIBD =  rbind(RIBD,data.frame(PT=p,RIBD_p))
  }
  RIBD_df = data.frame(PT   = rep(RIBD$PT,dim(RIBD)[2]-2),
                       Block_Name = rep(RIBD$Block,dim(RIBD)[2]-2),
                       GID = unlist(RIBD[,-(1:2)]))
  list(RIBD_F1= RIBD,RIBD_F2 =  RIBD_df)
}

PT_tst_RIBD_ls_f<-function(dat,GID='GID',Block_Name='Loc',nP=10,nBa=2,ss=1e4)
{
  set.seed(ss)
  RIBD = PT_Tr_RIBD_f(dat,GID,Block_Name,nP,nBa,ss)
  PT_tr =  RIBD$RIBD_F2
  Pos_tr =  match(paste(PT_tr$Block_Name,PT_tr$GID,sep='_'),
                  paste(dat[,Block_Name],dat[,GID],sep='_'))
  PT_tr$Pos_tr  =  Pos_tr
  PT_tst_ls =  vector('list',nP)
  for(p in 1:nP)
  {
    PT_tr_p =  PT_tr[PT_tr$PT==p,]
    PT_tst_ls[[p]] = setdiff(1:dim(dat)[1],PT_tr_p$Pos_tr)
  }
  names(PT_tst_ls) = paste0('PT',1:nP)
  PT_tst_ls
}


#==================================================================================
#Functions to be use to conform the partitions with the positions of testing data
#==================================================================================
#PT_tst_RIBD_ls_f(dat,GID='GID',Block_Name='Env',nP=10,nBa=2)
#and
#PT_tst_IBD_ls_f(dat,GID='GID',Block_Name='Env',nP=10,nBa=2)
#nP: No. of partitions to testing
#nBa: No of Blocks which each GID to be present in training data