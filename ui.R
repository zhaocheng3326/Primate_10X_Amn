#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram



shinyUI(
  navbarPage(
    "My Application",
    tabPanel("UMAP plot",
             fluidPage(
               # Application title
               titlePanel("UMAP and expression"),
               sidebarLayout(
                 sidebarPanel(
                   selectInput("textInput_days",label="select days(not not used currently)",choice=list('all'="All days","d10"="Day10","d12"="Day12","d14"="Day14"),selected="all"),
                   selectInput("textInput_ct",label="select type of expression",choice=list("SCT"="SCT","logTPM"="logTPM"),selected="1"),
                   selectizeInput("textInput_selG",label="select gene(GeneID_ENSID(like SOX2_ENSMFAG00000027781))",selected="SOX2_ENSMFAG00000027781",choice="SOX2_ENSMFAG00000027781"),
                   actionButton("actBshowExp", label ="Show exression")
                 ),
                 mainPanel(
                   fluidRow(
                     h3("umap plot",align="center")
                   ),
                   fluidRow(
                     column(3,plotOutput("PlotUMAP_full_d10")),
                     column(1,h1("|")),
                     column(3,plotOutput("PlotUMAP_full_d12")),
                     column(1,h1("|")),
                     column(3,plotOutput("PlotUMAP_full_d14"))
                   ),
                   fluidRow(
                     column(3,plotOutput("PlotUMAP_ICM_d10")),
                     column(1,h1("|")),
                     column(3,plotOutput("PlotUMAP_ICM_d12")),
                     column(1,h1("|")),
                     column(3,plotOutput("PlotUMAP_ICM_d14"))
                   ),
                   fluidRow(
                     h3("Feature plot",align="center")
                   ),
                   
                   fluidRow(
                     column(3,plotOutput("PlotFt_full_d10")),
                     column(1,h1("")),
                     column(3,plotOutput("PlotFt_full_d12")),
                     column(1,h1("")),
                     column(3,plotOutput("PlotFt_full_d14"))
                   ),
                   fluidRow(
                     column(3,plotOutput("PlotFt_ICM_d10")),
                     column(1,h1("")),
                     column(3,plotOutput("PlotFt_ICM_d12")),
                     column(1,h1("")),
                     column(3,plotOutput("PlotFt_ICM_d14"))
                   ),
                   fluidRow(
                     h3("violin plot",align="center")
                   ),
                   fluidRow(
                     column(3,plotOutput("PlotVln_full_d10")),
                     column(1,h4("")),
                     column(3,plotOutput("PlotVln_full_d12")),
                     column(1,h4("")),
                     column(3,plotOutput("PlotVln_full_d14"))
                   ),
                   fluidRow(
                     column(3,plotOutput("PlotVln_ICM_d10")),
                     column(1,h4("")),
                     column(3,plotOutput("PlotVln_ICM_d12")),
                     column(1,h4("")),
                     column(3,plotOutput("PlotVln_ICM_d14"))
                   )
               ),
               
               )
             )
    )
  )
)
