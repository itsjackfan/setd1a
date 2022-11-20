## RidgePlot and KS test on pseudotime
pseudo <- pseudotime(cds)
# add pseudotime to Seurat metadata
sdata <- AddMetaData(
  object = sdata,
  metadata = pseudo,
  col.name = 'pseudotime'
)
head(x = sdata[[]])

## define ridge plot function
drawRidge <- function(so){
  # bY CELL TYPE
  print(
    RidgePlot(
      so,
      features = "pseudotime",
      cols = NULL,
      idents = NULL,
      sort = FALSE,
      assay = NULL,
      group.by = 'nowakowski_med',
      y.max = 10,
      same.y.lims = TRUE,
      log = FALSE,
      ncol = NULL,
      slot = "data",
      stack = FALSE,
      combine = FALSE,
      fill.by = "ident"
    )
  )
  
  # By genotype
  print(
    RidgePlot(
      so,
      features = "pseudotime",
      cols = NULL,
      idents = NULL,
      sort = TRUE,
      assay = NULL,
      group.by = 'condition',
      #  y.max = 10,
      same.y.lims = TRUE,
      log = FALSE,
      ncol = NULL,
      slot = "data",
      stack = FALSE,
      combine = TRUE,
      fill.by = NULL
    ) + geom_density_ridges(scale = FALSE, alpha = 0.2)
  )
}

## extract sub cell populations if needed
sdata_sub = subset(x = sdata, subset = nowakowski_med %in% c("EN-V1","RG-early") )

pdf(file = "ridgeplots.pdf")
drawRidge(sdata)
drawRidge(sdata_sub)
dev.off()

## STATISTICS K-S test
ridgeData = sdata@meta.data

case=subset(ridgeData, ridgeData$condition =="FS")
control=subset(ridgeData, ridgeData$condition=="CT")

ridgeKS_res = ks.test(case$pseudotime , control$pseudotime)

write(unlist(ridgeKS_res),file="ridgeKS_res.txt")
# create ECDF of data
sample1 = control$pseudotime
sample2 = case$pseudotime
group <- c(rep("CT", length(sample1)), rep("FS", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1)
cdf2 <- ecdf(sample2)
# find min and max statistics to draw line between points of greatest distance
minMax <- seq(min(sample1, sample2), max(sample1, sample2), length.out=length(sample1))
x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )]
y0 <- cdf1(x0)
y1 <- cdf2(x0)

pdf(file = "Ks_ridgeplot.pdf")
ggplot(dat, aes(x = KSD, group = group, color = group))+
  stat_ecdf(size=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("pseudotime") +
  ylab("ECDF") +
  #geom_line(size=1) +
  geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
               linetype = "dashed", color = "red") +
  geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=8) +
  geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=8) +
  ggtitle("K-S Test: CT / FS") +
  theme(legend.title=element_blank())
dev.off()