mat = readRDS("Expression_affy_v2_data_matrix.rds")
zmat = t(scale(t(mat)))
immune = readRDS("immune_score_matrix_vZScoreAverages.rds")
### Unsupervised clustering
vars = apply(mat,1,sd)^2
# Threshold vector
lx = seq(3,50,1) / 10

pve = len = rep(NA,length(lx))

for(i in 1:length(lx)){
  #print(i)
  SDCUT = lx[i]
  len[i] = length(which(vars > SDCUT))
  
  zmat2 = zmat[which(vars > SDCUT),]
  pc = prcomp(t(zmat2))
  x3 = pc$x[,1:3]
  pve[i] = sum(apply(x3,2,sd)^2) / sum(apply(pc$x,2,sd)^2)
}#end for i

resp = as.factor(sapply(strsplit(colnames(zmat),"_"),"[",c(F,T)))
respCol = rep("red",length(resp))
respCol[which(resp=="LT")] = "turquoise4"
###################################################################
library(shiny)

shinyUI(
    tabsetPanel(id="Unsup",
                tabPanel(title="PCA animation",
                         fluidPage(
                           fluidRow(
                             column(12, align="left",
                                    sliderInput("animation", "Probe set variance cutoff:", 0.3, 5, 0.3,step = 0.1, animate=animationOptions(interval=500, loop=FALSE)))
                           ),
                           fluidRow(
                             column(6,plotOutput("vardensity",height = 350, width = 420)),
                             column(6,plotOutput("percentVar",height = 350, width = 420))
                           ),
                           fluidRow(
                             column(6,plotOutput("pc1pc2",height = 350, width = 420)),
                             column(6,plotOutput("pc1pc3",height = 350, width = 420))
                           )
                         )),
                tabPanel(title="Probe set table",
                         fluidPage(
                           fluidRow(
                             column(12,DT::dataTableOutput("geneTable"))
                           )
                         )
                ),
                tabPanel(title="Heat map",
                         fluidPage(
                           fluidRow(
                             column(4, align="left",uiOutput("slideHeatmap")),
                             column(5),
                             column(2, align="center",verbatimTextOutput("numGene"))
                           ),
                           fluidRow(
                             column(12,imageOutput("heatmap"))
                           ),
                           tags$style(type='text/css', "#numGene { width:100%; margin-top: 25px;}")
                         )
                )
    )
)
    