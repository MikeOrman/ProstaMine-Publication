oncoprint<-function(gene1,gene2,gene3,CNA){
  a<-CNA[c(gene1,gene2,gene3),]
  a<-as.data.frame(a)
  a[a == 0] <- ""
  a[a==-2] <- "Deep Deletion"
  a[a==-1] <- "Shallow Deletion"
  a[a==1] <- "Gain"
  a[a==2] <- "Amplification"
  a[4,]<-"XYZ"
  alter_fun <- list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
                gp = gpar(fill = "#CCCCCC", col = NA))
    },
    "XYZ"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="orange", col=NA))
    },
    "Shallow Deletion"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="light blue", col=NA))
    },
    "Deep Deletion"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="blue", col=NA))
    },
    "Gain"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="pink", col=NA))
    },
    "Amplification"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="red", col=NA))
    }
  )
  cols<-c("No alteration"="grey","Shallow Deletion"="light blue","Deep Deletion"="blue","Gain"="pink","Amplification"="red")
  ComplexHeatmap::oncoPrint(a[c(1:3),],alter_fun = alter_fun,col = cols,
                            row_names_gp = gpar(fontsize = 20), pct_gp = gpar(fontsize = 20),
                            heatmap_legend_param = list(title = "Alterations",
                                                        title_gp = gpar(fontsize=14),
                                                        grid_height = unit(8, "mm"),
                                                        grid_width = unit(8, "mm")))
}

oncoprint1<-function(gene1,gene3,CNA){
  a<-CNA[c(gene1,gene3),]
  a<-as.data.frame(a)
  a[a == 0] <- ""
  a[a==-2] <- "Deep Deletion"
  a[a==-1] <- "Shallow Deletion"
  a[a==1] <- "Gain"
  a[a==2] <- "Amplification"
  a[4,]<-"XYZ"
  alter_fun <- list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
                gp = gpar(fill = "#CCCCCC", col = NA))
    },
    "XYZ"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="orange", col=NA))
    },
    "Shallow Deletion"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="light blue", col=NA))
    },
    "Deep Deletion"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="blue", col=NA))
    },
    "Gain"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="pink", col=NA))
    },
    "Amplification"=function(x, y, w, h)
    {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), gp=gpar(fill="red", col=NA))
    }
  )
  cols<-c("No alteration"="grey","Shallow Deletion"="light blue","Deep Deletion"="blue","Gain"="pink","Amplification"="red")
  ComplexHeatmap::oncoPrint(a[c(1:2),],alter_fun = alter_fun,col = cols,
                            row_names_gp = gpar(fontsize = 20), pct_gp = gpar(fontsize = 20),
                            heatmap_legend_param = list(title = "Alterations",
                                                        title_gp = gpar(fontsize=14),
                                                        grid_height = unit(8, "mm"),
                                                        grid_width = unit(8, "mm")))
}