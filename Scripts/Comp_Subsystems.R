
a <- 0
for (i in subsys) {
  a <- a +1
  i <- subsys[1]
  reacts <- Info_React$ID[Info_React$Subsystem == i]
  name_sub <- paste0("Sub_", i)
  assign(name_sub, reacts)
}

## Import names of subsystems
Info_React <- read.delim("../reacts_subsystems.txt")
subsys <- unique(Info_React$subsystem)

## Import the WT FILE
mtfdir <- '../obj_PA14_Biomass/iPau21+allbidi.tab/'
name_file_mtf <- grep(".Rlog", list.files(mtfdir), value = TRUE)
file_wt_mtf <- paste(mtfdir, name_file_mtf, sep='/')

wt_mtf <- read.table(file_wt_mtf)
if (dim(wt_mtf)[2] == 3) {
  wt_mtf <- wt_mtf[,2:3]
}
colnames(wt_mtf) <- c('Reaction', 'Flux')

## Add name of subsystem to reactions and order the reactions by its subsystem
wt_mtf$Subsystem <- Info_React$Subsystem
wt_mtf_Sub <- wt_mtf[order(wt_mtf$Subsystem),]

## Obtain the reactions of each subsystem
subsys <- levels(wt_mtf_Sub$Subsystem)
subsys[1]

wt_mtf_Sub[levels(wt_mtf_Sub$Subsystem)[1],]

reacts <- Info_React$ID[Info_React$Subsystem == subsys[1]]
