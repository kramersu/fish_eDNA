library(dada2)
library(ShortRead)
library(Biostrings)
library(reticulate)
library(dplyr)
library(phyloseq)

setwd("C:/Users/KRAMERS/Documents/Gerald_fish_eDNA/raw/12S/")
#rename files
#myfiles<-list.files(pattern="fastq.gz")
#gsub("23-002232-","",myfiles)
#file.rename(from=myfiles,to=gsub("23-002232-","",myfiles))

myfiles<-list.files(pattern="fastq.gz")

path_raw<-c("C:/Users/KRAMERS/Documents/Gerald_fish_eDNA/raw/12S/")
path<-c("C:/Users/KRAMERS/Documents/Gerald_fish_eDNA/12S/")
fnFs<-sort(list.files(path_raw,pattern="_R1_001.fastq.gz",full.names=TRUE))
fnRs<-sort(list.files(path_raw,pattern="_R2_001.fastq.gz",full.names=TRUE))

#Find the primers in the raw data
FWD<-"CCGGTAAAACTCGTGCCAGC"
REV<-"CATAGTGGGGTATCTAATCCCAGTTTG"


allOrients<-function(primer){
  require(Biostrings)
  dna<-DNAString(primer)
  orients<-c(Forward=dna,Complement=complement(dna),Reverse=reverse(dna),RevComp=reverseComplement(dna))
  return(sapply(orients,toString))
}

FWD.orients<-allOrients(FWD)
REV.orients<-allOrients(REV)



primerHits<-function(primer,fn){
  nhits<-vcountPattern(primer,sread(readFastq(fn)),fixed=FALSE)
  return(sum(nhits>0))
}

rbind(FWD.ForwardReads=sapply(FWD.orients,primerHits,fnFs[[1]]),
      FWD.ReverseReads=sapply(FWD.orients,primerHits,fnRs[[1]]),
      REV.ForwardReads=sapply(REV.orients,primerHits,fnFs[[1]]),
      REV.ReverseReads=sapply(REV.orients,primerHits,fnRs[[1]])
)

cutadapt<-"C:/Users/KRAMERS/AppData/Roaming/Python/Python311/Scripts/cutadapt"
system2(cutadapt,args='--version')

path.cut<-file.path(path,"cutadapt")
if(!dir.exists(path.cut))dir.create(path.cut)
fnFs.cut<-file.path(path.cut,basename(fnFs))
fnRs.cut<-file.path(path.cut,basename(fnRs))

FWD.RC<-dada2:::rc(FWD)
REV.RC<-dada2:::rc(REV)

R1.flags<-paste("-g",FWD,"-a",REV.RC)
R2.flags<-paste("-G",REV,"-A",FWD.RC)

#Trim primers with cutadapt
for(i in seq_along(fnFs)){
  system2(cutadapt, arg =c(R1.flags,R2.flags,"-n",2,"-m",10,"-o",fnFs.cut[i],"-p",fnRs.cut[i],fnFs[i],fnRs[i]))
}

rbind(FWD.ForwardReads=sapply(FWD.orients,primerHits,fnFs.cut[[1]]),
      FWD.ReverseReads=sapply(FWD.orients,primerHits,fnRs.cut[[1]]),
      REV.ForwardReads=sapply(REV.orients,primerHits,fnFs.cut[[1]]),
      REV.ReverseReads=sapply(REV.orients,primerHits,fnRs.cut[[1]])
)

cutFs<-sort(list.files(path.cut,pattern="_R1_001.fastq.gz",full.names=TRUE))
cutRs<-sort(list.files(path.cut,pattern="_R2_001.fastq.gz",full.names=TRUE))
get.sample.name<-function(fname)strsplit(basename(fname),"_L001")[[1]][1]
sample.names<-unname(sapply(cutFs,get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[1:6])#read quality drops after about 225bp
plotQualityProfile(cutRs[1:6])#read quality drops around 225bp

filtFs<-file.path(path,"filtered",basename(cutFs))
filtRs<-file.path(path,"filtered",basename(cutRs))

out<-filterAndTrim(cutFs,filtFs,cutRs,filtRs,maxN=0,maxEE=c(2,2),truncQ = 2,
                   minLen=50,rm.phix=TRUE,compress=TRUE)

head(out)
exists<-file.exists(filtFs)#all files have reads

errF<-learnErrors(filtFs[exists],multithread = TRUE)
errR<-learnErrors(filtRs[exists],multithread = TRUE)

plotErrors(errF,nominalQ=TRUE)
plotErrors(errR,nominalQ=TRUE)

derepFs<-derepFastq(filtFs[exists],verbose=TRUE)
derepRs<-derepFastq(filtRs[exists],verbose=TRUE)

dadaFs<-dada(derepFs,err=errF,multithread=TRUE,pool="pseudo")
dadaRs<-dada(derepRs,err=errR,multithread=TRUE,pool="pseudo")

mergers<-mergePairs(dadaFs,derepFs,dadaRs,derepRs,verbose=TRUE)
names(mergers)<-sample.names[exists]
seqtab<-makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))#most sequences are 253 bp long
seqtab2<-seqtab[,nchar(colnames(seqtab)) %in% 250:281]

seqtab2.nochim<-removeBimeraDenovo(seqtab2,method="consensus",multithread=TRUE,verbose=TRUE)
sum(seqtab2.nochim)/sum(seqtab2)#99.9% of reads remain

getN<-function(x) sum(getUniques(x))
out_ex<-row.names(out)[exists]
track<-cbind(out,sapply(dadaFs,getN),sapply(dadaRs,getN),sapply(mergers,getN),
             rowSums(seqtab2.nochim))
colnames(track)<-c("input","filtered","denoisedF","denoisedR","merged","nochim")

ps_base<-phyloseq(otu_table(seqtab2.nochim,taxa_are_rows = FALSE))
dna<-Biostrings::DNAStringSet(taxa_names(ps_base))
names(dna)<-taxa_names(ps_base)
ps_base<-merge_phyloseq(ps_base,dna)
ps_sub<-prune_taxa(taxa_sums(ps_base)>2,ps_base)#2958 taxa remain



taxa_names(ps_sub)<-paste0("ASV",seq(ntaxa(ps_sub)))
ps_sub %>%
  refseq() %>%
  Biostrings::writeXStringSet("C:/Users/KRAMERS/Documents/Gerald_fish_eDNA/12S/ASV.fna",append=FALSE,
                              compress=FALSE,compression_level=NA,format="fasta")


