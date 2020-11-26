#!/usr/bin/env Rscript
#Autor: Kenny Pinheiro
#Data: 26/04/2020

library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Multifasta or fasta file name", metavar="character"),
  make_option(c("-n", "--nuc"), type = "integer", default = 3,
              help = "Number of nucleotides in repetitions [default %default]",
              metavar = "number"),
  make_option(c("-l", "--labels"), type = "integer", default = 3,
              help = "X-Axis-Labels - [1] will show all labels. Higher values will show less [default %default]",
              metavar = "number"),
  make_option(c("-f", "--font"), type = "integer", default = 10,
              help = "Font Size of X-Axis [default %default]",
              metavar = "number"),
  make_option(c("-o", "--out"), type="character", default = "out.png", 
              help="output file name [default = %default]", metavar="character")
); 
  
opt_parser = OptionParser(usage = "%prog -i file.fasta [options]"
                          ,option_list=option_list);
opt = parse_args(opt_parser);
  
if(is.null(opt$out)) opt$out <- gsub(".fasta$", ".png", opt$input)


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input fasta file)", call.=FALSE)
}

library(Biostrings)
library(ggplot2)

dna = readDNAStringSet(opt$input)
nomes <- names(dna)
tam <- length(dna)
df <- List()

for(i in 1:tam){
  df[[i]] <- as.data.frame(oligonucleotideFrequency(dna[i],opt$nuc))
  Repetitions <- names(df[[i]])
  df[[i]] <- as.data.frame(t(df[[i]]))
  df[[i]][,"Repetitions"] <- Repetitions
  df[[i]][,"Genomes"] <- rep(nomes[i],length(df[[i]]$V1))
  names(df[[i]]) <- c("Freq","Repetitions","Genomes")
}

DF <- as.data.frame(df)
DF$group_name <- NULL

##Modify width and height values to scale the output PNG file
png(filename=opt$out, width= 1024, height= 800)

ggplot(DF, aes(x=Repetitions, y=Freq, group=Genomes, fill=Genomes)) + 
  geom_line(aes(color=Genomes)) + 
  geom_point(aes(color=Genomes, shape=Genomes)) + 
  theme(legend.direction='vertical', legend.position="right") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size=opt$font)) +

  scale_x_discrete(
    labels = Repetitions[seq(1,length(Repetitions),opt$labels)],
    breaks = Repetitions[seq(1,length(Repetitions),opt$labels)]
  ) 

dev.off()


#If the user wanna got tables of frequency for each sample, please remove the '#' from the lines below:

#library(dplyr)
#data <- list()
#for (i in 1:tam){
#  data[[i]] <- DF %>% filter(Genomes == names((table(as.factor(DF$Genomes))[i])))
#  write.table(data[[i]][,2:3], file=data[[i]][1,4], sep="\t", row.names = FALSE, quote = FALSE)
#}