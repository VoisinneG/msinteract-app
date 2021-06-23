bait_list <- c("Grb2", "Cbl", "Cblb", "Fyb", "Inpp5d", "Itk", "Lck", "Lcp2", "Nfatc2",
               "Ptpn22", "Vav1", "Plcg1", "Themis", "Ptpn6", "Nck1")

name_dir <- "~/ownCloud/++Work/++Research/++Projects/Integration_15_interactomes/Analysis/InteRactomes/"

for(bait in bait_list){
  
  load(paste(name_dir,"Interactome_", bait,".rda",sep=""))
  Interactome <- get(paste("Interactome_", bait,sep=""))
  assign(paste("Interactome_", bait, sep=""), Interactome)
  save(list=paste("Interactome_", bait, sep=""), file=paste("./data/Interactome_", bait,".rda", sep=""))
  
}


#Load proteomes
load("~/ownCloud/++Work/++Research/++Projects/Proteomes/Merge_Proteome_OST_and_CD4_kinetics/output/proteome_CD4.rda")
save(proteome_CD4, file = "./data/proteome_CD4.rda")

load("~/ownCloud/++Work/++Research/++Projects/Proteomes/Proteome_CD4+_transfected_expanded/output/proteome_CD4_expanded.rda")
save(proteome_CD4_expanded, file = "./data/proteome_CD4_expanded.rda")

load("~/ownCloud/++Work/++Research/++Projects/Proteomes/Proteome_Jurkat_Itsn2_KO/output/proteome_Jurkat.rda")
save(proteome_Jurkat, file = "./data/proteome_Jurkat.rda")
