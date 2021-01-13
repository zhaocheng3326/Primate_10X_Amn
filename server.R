#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
rm(list = ls())
suppressMessages(library(shiny))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
#suppressMessages(library(rgl))
#suppressMessages(library(shinysky))
suppressMessages(library(plotly))

source("quick.functions.R")

if (file.exists("data/process.data.Rdata")) {
    #load("data/process.data.Rdata",verbose = T)    ### full dataset
  load("data/process.smalldata.Rdata",verbose = T)  ### smaller one for test
}else{
    process.data <- FunPreProcessData(TRUE)

    ENS.ID <- process.data$GI
    GIEN <- data.frame(ENS=ENS.ID,GI=Funid_translate(ENS.ID)) %>% tbl_df() %>% mutate(ENS=as.vector(ENS),GI=as.vector(GI)) %>% mutate(SGI=ifelse(GI==ENS,as.vector(ENS),as.vector(paste(GI,ENS,sep="_"))))

    data.out <- process.data$data
    meta.out <- process.data$meta

   save(data.out,meta.out,ENS.ID,GIEN,file="data/process.data.Rdata")
   #save(data.out,meta.out,ENS.ID,GIEN,file="data/process.smalldata.Rdata")
}




#' loading data
# data=data.out
# meta=meta.out
# day="d10"
# type="ICM"
# ct="logTPM"
# GIEN=GIEN
# sid="ASB6_ENSMFAG00000007597"
# exp.limit=3.5
# scale=TRUE
# gene.sel.scal.exp <- FunQuickGetData(data=data.out,meta=meta.out,day="d10",type="ICM",ct="logTPM",gien=GIEN,sid="ASB6_ENSMFAG00000007597",exp.limit=3.5,scale=FALSE)
# 

temp.plot <- list()
temp.plot$PlotUMAP_ICM_d10 <- FunUMAPPlot(meta.out$ICM$d10,paste("ICM","d10"))
temp.plot$PlotUMAP_full_d10 <- FunUMAPPlot(meta.out$full$d10,paste("All data","d10"))
temp.plot$PlotUMAP_ICM_d12 <- FunUMAPPlot(meta.out$ICM$d12,paste("ICM","d12"))
temp.plot$PlotUMAP_full_d12 <- FunUMAPPlot(meta.out$full$d12,paste("All data","d12"))
temp.plot$PlotUMAP_ICM_d14 <- FunUMAPPlot(meta.out$ICM$d14,paste("ICM","d14"))
temp.plot$PlotUMAP_full_d14 <- FunUMAPPlot(meta.out$full$d14,paste("All data","d14a"))

shinyServer(function(input, output, session) {
  #day="d10"
  #type="ICM"
  #sid="PGBD2_ENSMFAG00000044637"
  updateSelectizeInput(session,'textInput_selG', choices = GIEN$SGI, server = TRUE)
  for ( n in names(temp.plot)) {
    print(n)
    output[[n]] <- renderPlot(temp.plot[[n]])
  }

 
  output$PlotUMAP_full_d10 <- renderPlot(temp.plot$PlotUMAP_full_d10)
  output$PlotUMAP_ICM_d10 <- renderPlot(temp.plot$PlotUMAP_ICM_d10)
  output$PlotUMAP_full_d12 <- renderPlot(temp.plot$PlotUMAP_full_d12)
  output$PlotUMAP_ICM_d12 <- renderPlot(temp.plot$PlotUMAP_ICM_d12)
  output$PlotUMAP_full_d14 <- renderPlot(temp.plot$PlotUMAP_full_d14)
  output$PlotUMAP_ICM_d14 <- renderPlot(temp.plot$PlotUMAP_ICM_d14)
  
  
  plotGeneExpevent <- eventReactive(input$actBshowExp, {
    para <- list(ct=input$textInput_ct,sid=as.character(input$textInput_selG))
    para
  })
    
  
    output$PlotFt_ICM_d10 <- renderPlot({para <-plotGeneExpevent();FunSelFeaturePlot(data=data.out,meta=meta.out,day="d10",type="ICM",ct=para$ct,gien=GIEN,sid=para$sid,plot.col=c("grey","blue"))})
    output$PlotFt_full_d10 <- renderPlot({para <-plotGeneExpevent();FunSelFeaturePlot(data=data.out,meta=meta.out,day="d10",type="full",ct=para$ct,gien=GIEN,sid=para$sid,plot.col=c("grey","blue"))})
    output$PlotFt_ICM_d12 <- renderPlot({para <-plotGeneExpevent();FunSelFeaturePlot(data=data.out,meta=meta.out,day="d12",type="ICM",ct=para$ct,gien=GIEN,sid=para$sid,plot.col=c("grey","blue"))})
    output$PlotFt_full_d12 <- renderPlot({para <-plotGeneExpevent();FunSelFeaturePlot(data=data.out,meta=meta.out,day="d12",type="full",ct=para$ct,gien=GIEN,sid=para$sid,plot.col=c("grey","blue"))})
    output$PlotFt_ICM_d14 <- renderPlot({para <-plotGeneExpevent();FunSelFeaturePlot(data=data.out,meta=meta.out,day="d14",type="ICM",ct=para$ct,gien=GIEN,sid=para$sid,plot.col=c("grey","blue"))})
    output$PlotFt_full_d14 <- renderPlot({para <-plotGeneExpevent();FunSelFeaturePlot(data=data.out,meta=meta.out,day="d14",type="full",ct=para$ct,gien=GIEN,sid=para$sid,plot.col=c("grey","blue"))})
   
    
    output$PlotVln_ICM_d10 <- renderPlot({para <-plotGeneExpevent();FunSelVlnePlot(data=data.out,meta=meta.out,day="d10",type="ICM",ct=para$ct,gien=GIEN,sid=para$sid)})
    output$PlotVln_full_d10 <- renderPlot({para <-plotGeneExpevent();FunSelVlnePlot(data=data.out,meta=meta.out,day="d10",type="full",ct=para$ct,gien=GIEN,sid=para$sid)})
    output$PlotVln_ICM_d12 <- renderPlot({para <-plotGeneExpevent();FunSelVlnePlot(data=data.out,meta=meta.out,day="d12",type="ICM",ct=para$ct,gien=GIEN,sid=para$sid)})
    output$PlotVln_full_d12 <- renderPlot({para <-plotGeneExpevent();FunSelVlnePlot(data=data.out,meta=meta.out,day="d12",type="full",ct=para$ct,gien=GIEN,sid=para$sid)})
    output$PlotVln_ICM_d14 <- renderPlot({para <-plotGeneExpevent();FunSelVlnePlot(data=data.out,meta=meta.out,day="d14",type="ICM",ct=para$ct,gien=GIEN,sid=para$sid)})
    output$PlotVln_full_d14 <- renderPlot({para <-plotGeneExpevent();FunSelVlnePlot(data=data.out,meta=meta.out,day="d14",type="full",ct=para$ct,gien=GIEN,sid=para$sid)})
    
    #output$PlotVln_ICM_d10 <- renderPlot(FunSelVlnePlot(data=data.out,meta=meta.out,day="d10",type="ICM",ct=para$ct,gien=GIEN,sid=sid))
    #output$PlotVln_full_d10 <- renderPlot(FunSelVlnePlot(data=data.out,meta=meta.out,day="d10",type="full",ct=para$ct,gien=GIEN,sid=sid))
    #output$PlotVln_ICM_d12 <- renderPlot(FunSelVlnePlot(data=data.out,meta=meta.out,day="d12",type="ICM",ct=para$ct,gien=GIEN,sid=sid))
    #output$PlotVln_full_d12 <- renderPlot(FunSelVlnePlot(data=data.out,meta=meta.out,day="d12",type="full",ct=para$ct,gien=GIEN,sid=sid))
    #output$PlotVln_ICM_d14 <-  renderPlot(FunSelVlnePlot(data=data.out,meta=meta.out,day="d14",type="ICM",ct=para$ct,gien=GIEN,sid=sid))
    #output$PlotVln_full_d14 <- renderPlot(FunSelVlnePlot(data=data.out,meta=meta.out,day="d14",type="full",ct=para$ct,gien=GIEN,sid=sid))
    

}
)
                    
