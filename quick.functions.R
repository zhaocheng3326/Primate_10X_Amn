# script to show the plot

rename <- dplyr::rename
filter <- dplyr::filter

FunPreProcessData <- function(OV=FALSE) {
  if (OV) {
    data.out <- list()
    meta.out <- list()
    GI <- c()
    for (type in c("ICM","full")) {
      for (day in c("d10","d12","d14")) {
        for (ct in c("raw","SCT")) {
          temp  <- read.delim(paste0("raw_data/",type,"_dataset/WT_",day,".data.",type,".",ct,".txt"),stringsAsFactors = F,head=T,row.names = 1) 
          GI <- c(GI,rownames(temp))
          if (ct =="raw") {
            temp <- apply(temp,2,function(x){x/sum(x)*1000000})%>% log1p() %>% as.data.frame()
            data.out[[type]][[day]][["logTPM"]] <- temp %>% tibble::rownames_to_column("ENS") %>% tbl_df() 

          }else{
            data.out[[type]][[day]][[ct]] <- temp %>% tibble::rownames_to_column("ENS") %>% tbl_df()
          }
        }
        temp.meta <- data.frame("cell"=colnames(temp),lineage=read.delim(paste0("raw_data/",type,"_dataset/WT_",day,".data.",type,".metadata.txt"),stringsAsFactors = F,head=F)$V2) %>% tbl_df() %>% inner_join(read.delim(paste0("raw_data/",type,"_dataset/WT_",day,".data.",type,".umap.txt"),stringsAsFactors = F,head=T,row.names = 1)  %>% tibble::rownames_to_column("cell") %>% mutate(cell=paste0("X",cell)),by="cell" ) %>% mutate(lineage=as.vector(lineage))
        meta.out[[type]][[day]] <- temp.meta
      }
    }
  }
  return(list(data=data.out,meta=meta.out,GI=unique(GI)))
}

Funid_translate <- function(id) {
  annotations <- read.table('./raw_data/mart_export_idtranslate_fun.txt',header = T, sep = ",", stringsAsFactors = F)
  translated <- vector(mode = 'character',length = length(id)) 
  for(i in 1:length(id)){
    if(!is.na(match(id[i], annotations[,1]))){
      if(annotations[match(id[i], annotations[,1]),2] != ''){
        translated[i] <- annotations[match(id[i], annotations[,1]),2]
      }
      else{
        translated[i] <- id[i]
      }
      if(length(grep('LOC', translated[i])!=0)){
        if(annotations[match(id[i], annotations[,1]),3] != ''){
          translated[i] <- annotations[match(id[i], annotations[,1]),3]
        }
      }
    }
    else if (!is.na(match(id[i], annotations[,2]))){
      translated[i] <- annotations[match(id[i], annotations[,2]),1]
    }
    else if (!is.na(match(id[i], annotations[,3]))){
      translated[i] <- annotations[match(id[i], annotations[,3]),1]
    }
    else{
      translated[i] <- id[i]
      warning(paste0('Gene ',id[i], ' not found'))
    }
  }
  return(translated)
}
FunTitle <- function() {
  p <- theme(plot.title = element_text(hjust=0.5,face="bold"))
  return(p)
}


FunEmptyPlot <- function(GI, DS) {
  return(
    ggplot(data.frame(1,1))+geom_blank()+theme_void()+ggtitle(paste(GI,"isn't in", DS,sep="\n"))+theme(plot.title = element_text(hjust=0.5,face="bold"))
  )
}

FunColSet <- function() {
  p <- scale_color_manual(values=c("Amn-Progenitors"="#FCCDE5","Amnion"="#FCCDE5","Amnion 1"="#FCCDE5","Amnion 2"="#BC80BD","EPI"="#a84c28","EPI + EPI-derived"="#E41A1C","ExE-Mesenchyme"="#377EB8","ExE-Mesoderm"="#377EB8","Transition"="#FD6A02","Trophoblast"="#984EA3","VE/YE"="#4DAF4A","emerging ExE-Mesoderm"="#BEBADA","emerging mesoderm 1"="#70C1B3","emerging mesoderm 2"="#7a7a79"))
  return(p)
}
FunFillColSet <- function() {
  p <- scale_fill_manual(values=c("Amn-Progenitors"="#FCCDE5","Amnion"="#FCCDE5","Amnion 1"="#FCCDE5","Amnion 2"="#BC80BD","EPI"="#a84c28","EPI + EPI-derived"="#E41A1C","ExE-Mesenchyme"="#377EB8","ExE-Mesoderm"="#377EB8","Transition"="#FD6A02","Trophoblast"="#984EA3","VE/YE"="#4DAF4A","emerging ExE-Mesoderm"="#BEBADA","emerging mesoderm 1"="#70C1B3","emerging mesoderm 2"="#7a7a79"))
  return(p)
}

FunUMAPPlot <- function(meta,DS) {
  meta <- meta
  temp.text.pos <-  meta %>%group_by(lineage) %>% summarise(x=median(UMAP_1),y=median(UMAP_2)) %>% ungroup()
  temp.plot <-  ggplot(meta)+geom_point(mapping=aes(x=UMAP_1,y=UMAP_2,col=lineage))+theme_void()+FunColSet()+geom_text(temp.text.pos,mapping=aes(x=x,y=y,label=lineage))+ggtitle(DS)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ theme(legend.position = "none")
  return(temp.plot)
}

FunQuickGetData <- function(data=data.out,meta=meta.out,day="d10",type="ICM",ct="logTPM",gien=GIEN,sid="ASB6_ENSMFAG00000007597",exp.limit=3.5,scale=TRUE) {
  ens.id <- gien %>% filter(SGI==sid) %>% pull(ENS)
  temp.exp <- data[[type]][[day]][[ct]]%>% filter(ENS==ens.id) %>% tibble::column_to_rownames("ENS")
  if (nrow(temp.exp)==0) {
    return(temp.exp[,1:2])
  }else{
    if (!scale) {
      temp.sel.exp <- t(apply(temp.exp,1,scale))
    }else{
      temp.sel.exp <- temp.exp
    }
    
    colnames(temp.sel.exp) <- colnames(temp.exp)
    rownames(temp.sel.exp) <- rownames(temp.exp)
    temp.sel.exp[temp.sel.exp>exp.limit] <- exp.limit
    temp.sel.exp[temp.sel.exp<  exp.limit*(-1)] <- -exp.limit*(-1)
    
    return(temp.sel.exp %>% tbl_df() %>% gather(cell,sv) %>%mutate(SGI=sid) %>% full_join(meta[[type]][[day]],by="cell"))
  }
  
}


FunSelFeaturePlot <- function (data=data.out,meta=meta.out,day="d10",type="ICM",ct="logTPM",gien=GIEN,sid="ASB6_ENSMFAG00000007597",plot.col=c("grey","blue"))  {
  gene.sel.scal <- FunQuickGetData (data=data,meta=meta,day=day,type=type,ct=ct,gien=gien,sid=sid,exp.limit=3.5,scale=TRUE)
  if (type=="full") {type="all dataset"}
  if (nrow(gene.sel.scal) ==0) {
    temp.plot <- FunEmptyPlot(sid,paste(day,type))
  }else{
    temp.plot <- gene.sel.scal %>% ggplot()+geom_point(mapping = aes(x=UMAP_1,y=UMAP_2,color=sv),size=1)+theme_void()+ggtitle(paste(sid,paste("expressed",day,type),sep="\n"))+theme(plot.title = element_text(hjust=0.5,face="bold"))+scale_color_gradientn(colors=plot.col)+theme(legend.position = "right",legend.title = element_blank())
    
  }
  return(temp.plot)
 
}
FunSelVlnePlot <- function (data=data.out,meta=meta.out,day="d10",type="ICM",ct="logTPM",gien=GIEN,sid="ASB6_ENSMFAG00000007597")  {
  gene.sel.scal <- FunQuickGetData (data=data,meta=meta,day=day,type=type,ct=ct,gien=gien,sid=sid,exp.limit=3.5,scale=TRUE)
  if (type=="full") {type="all dataset"}
  if (nrow(gene.sel.scal) ==0) {
    temp.plot <- FunEmptyPlot(sid,paste(day,type))
  }else{
    temp.plot <- gene.sel.scal %>% ggplot()+geom_violin(mapping = aes(x=lineage,y=sv,fill=lineage))+theme_classic()+ggtitle(paste(sid,paste("in",day,type),sep="\n"))+theme(plot.title = element_text(hjust=0.5,face="bold"))+FunFillColSet()+theme(legend.position = "none")+ylab(ct)+xlab("")+theme(axis.text.x=element_text(angle = 45)) 
  }
  return(temp.plot)
}
  
  
#data.out$ICM %>% lapply(function(x){dim(x$SCT)})
#meta.out$ICM %>% lapply(function(x){dim(x)})
FunPreProcessData_small <- function(OV=FALSE) {
  if (OV) {
    data.out <- list()
    meta.out <- list()
    GI <- c()
    for (type in c("ICM","full")) {
      for (day in c("d10","d12","d14")) {
        for (ct in c("raw","SCT")) {
          temp  <- read.delim(paste0("raw_data/",type,"_dataset/WT_",day,".data.",type,".",ct,".txt"),stringsAsFactors = F,head=T,row.names = 1) 
          GI <- c(GI,rownames(temp))
          if (ct =="raw") {
            temp <- apply(temp,2,function(x){x/sum(x)*1000000})%>% log1p() %>% as.data.frame()
            data.out[[type]][[day]][["logTPM"]] <- temp %>% tibble::rownames_to_column("ENS") %>% tbl_df() %>% head(3000)
            
          }else{
            data.out[[type]][[day]][[ct]] <- temp %>% tibble::rownames_to_column("ENS") %>% tbl_df()%>% head(3000)
          }
        }
        temp.meta <- data.frame("cell"=colnames(temp),lineage=read.delim(paste0("raw_data/",type,"_dataset/WT_",day,".data.",type,".metadata.txt"),stringsAsFactors = F,head=F)$V2) %>% tbl_df() %>% inner_join(read.delim(paste0("raw_data/",type,"_dataset/WT_",day,".data.",type,".umap.txt"),stringsAsFactors = F,head=T,row.names = 1)  %>% tibble::rownames_to_column("cell") %>% mutate(cell=paste0("X",cell)),by="cell" ) %>% mutate(lineage=as.vector(lineage))
        meta.out[[type]][[day]] <- temp.meta
      }
    }
  }
  return(list(data=data.out,meta=meta.out,GI=unique(GI)))
}

