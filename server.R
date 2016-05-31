library(shiny)

shinyServer(function(input, output) { 
  
  varcut <- reactive({
    input$animation
  })
  
  output$vardensity <- renderPlot({
    vc = varcut()
    dx = density(vars)
    
    par(mar=c(4,4,3,1)+0.1)
    plot(dx,log="x",lwd=3,col=4,main="Distribution of probe set variances",xlab="Probe set variance across 30 samples",cex.lab=1.5,cex.axis=1.5)
    abline(v=vc,col=2,lwd=3)
    mtext(text="Probe set var cutoff",side=3,line=0,at=vc,cex=1,col=2)
    
    ty = floor(max(dx$y))-0.5
    w = which(lx==vc)
    text(vc,ty,paste(len[w],"\nprobe sets",sep=""),pos=4,offset=0.3,col=2,font=2)
  })
  
  output$percentVar <- renderPlot({
    vc = varcut()
    w = which(lx==vc)
    par(mar=c(4,4,3,1)+0.1)
    plot(lx[1:w],100*pve[1:w],pch=16,type="p",xlab="Probe set variance cutoff",ylab="% variance explained",cex.lab=1.5,col=4,
         main="Percent variance explained by first 3 PCs",xlim=range(lx),ylim=c(100*min(pve),100*max(pve)),cex.axis=1.5,cex=1.5)
  })
  
  getScoreMat <- reactive({
    vc = varcut()
    zmat2 = zmat[which(vars > vc),]
    pc = prcomp(t(zmat2))
    pc$x
  })
  
  output$pc1pc2 <- renderPlot({
    x = getScoreMat()
    par(mar=c(4,4,3,1)+0.1)
    plot(x[,1],x[,2],pch=16,xlab="PC1",ylab="PC2",cex.lab=1.5,col=respCol,
         main="PC1 vs PC2",cex.axis=1.5,cex=1.5)
  })
  
  output$pc1pc3 <- renderPlot({
    x = getScoreMat()
    par(mar=c(4,4,3,1)+0.1)
    plot(x[,1],x[,3],pch=16,xlab="PC1",ylab="PC3",cex.lab=1.5,col=respCol,
         main="PC1 vs PC3",cex.axis=1.5,cex=1.5)
  })
  
  getTable <- reactive({
    vc = varcut()
    w = which(vars > vc)
    # significance between LT and ST
    pw = hiLT = rep(NA,length(w))
    for(i in 1:length(w)){
      pw[i] = wilcox.test(mat[w[i],] ~ resp,exact=FALSE)$p.value
      hiLT[i] = median(mat[w[i],which(resp=="LT")]) > median(mat[w[i],which(resp=="ST")])
    }#end for i
    
    # variance
    tmp = vars[w]
    # clean the XXX
    names(tmp) = gsub("XXX","",names(tmp))
    # data frame
    df = data.frame(PROBE.SET=names(tmp),VARIANCE=as.numeric(tmp),HIGH.in.LT=hiLT,PVAL=signif(pw,3))
    # order according to variance
    ordx = order(df[,2],decreasing=T)
    df = df[ordx,]
    uniqGene = sapply(strsplit(as.character(df[,1]),"//"),"[",c(T,F))
    url = paste0("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",uniqGene)
    df[,1] = paste0("<a href='",url,"' target='_blank'>",df[,1],"</a>")
    df
  })
  
  output$geneTable <- DT::renderDataTable({
    DT::datatable(getTable(),
                  escape=FALSE,
                  rownames = FALSE,
                  filter="top",
                  options = list(pageLength = 20),
                  caption = 'Probe sets sorted from highest to lowest variance')
  })
  
  output$slideHeatmap <- renderUI({
    vc = varcut()
    sliderInput("slideHeatmap", "Probe set variance cutoff:", 0.3, 5, vc,step = 0.1)
  })  
  
  output$heatmap <- renderImage({
    vc = input$slideHeatmap
    filename <- normalizePath(file.path('./heatmaps',paste0('cutoff',vc, '.png')))
    #filename <- normalizePath(file.path(paste0('cutoff',vc, '.png')))
    # Return a list containing the filename and alt text
    list(src = filename,
         width = 900,
         alt = paste("Cutoff", vc))
  }, deleteFile = FALSE)
  
  output$numGene <- renderText({
    w = which(lx==input$slideHeatmap)
    paste0(len[w],"\nprobe sets")
  })
}
)




