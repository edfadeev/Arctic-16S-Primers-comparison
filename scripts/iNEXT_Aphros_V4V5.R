library("iNEXT"); packageVersion("iNEXT")
library("ggplot2"); packageVersion("ggplot2")

PS85_BAC <- read.csv(file ="/scratch2/efadeev/Primers_comp/iNEXT/Data/V4V5_otu_table.csv", header = T)
rownames(PS85_BAC) <- PS85_BAC$X
PS85_BAC <- PS85_BAC[,-1]

iNEXT.out <- iNEXT(as.data.frame(PS85_BAC), q=1,
                   datatype="abundance", conf = 0.95, nboot = 50)

rare <-fortify(iNEXT.out, type=1)
write.csv(rare, file ="/scratch2/efadeev/Primers_comp/iNEXT/out/V4V5-1.csv")

rare <-fortify(iNEXT.out, type=2)
write.csv(rare, file ="/scratch2/efadeev/Primers_comp/iNEXT/out/V4V5-2.csv")

rare <-fortify(iNEXT.out, type=3)
write.csv(rare, file ="/scratch2/efadeev/Primers_comp/iNEXT/out/V4V5-3.csv")