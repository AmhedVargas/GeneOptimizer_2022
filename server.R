########GeneBuilder####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###Server

#Load libraries
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(Biostrings)
library(Cairo)
library(stringdist)
library(shinyjs)
library(DT)

##Multiple user related
library(promises)
library(future)
plan(multisession) 
options(future.rng.onMisuse="ignore")

#####Data##########################
##Codon Adapt
CAIS=read.table("DATA/Ultimate_aminos2.txt",sep="\t",header=T)
rownames(CAIS)= toupper(as.character(CAIS$Codon))
codons=unique(CAIS$Amino)

AAtoCodF=list()

for(i in 1:length(codons)){
  AAtoCodF=append(AAtoCodF, list(CAIS[which(CAIS$Amino == codons[i]),]))
  names(AAtoCodF)[i]=as.character(codons[i])
}

IntronSeqs=read.table("DATA/Introns.csv",sep=",",header=F,row.names=1)

##piRNAs
###There was no need to have them separated, If I make a RDs, the best will be to remove from the list these files
Pies=readLines("DATA/HengPies.txt")

PiesNA=readLines("DATA/HengNames.txt")

PiesFin=cbind(Pies,PiesNA)
rownames(PiesFin)=as.character(Pies)

##Enzymes
enzy=read.table("DATA/Enzymes.txt", sep="\t", colClasses = "character",header=T)
rownames(enzy)=as.character(enzy$Enzyme)
GoldenE=enzy[which(enzy$GG == "TRUE"),]
OtherE=enzy[-c(which(enzy$GG == "TRUE")),]

##TO DO: 
#It will be better to make the data transformations and save them as a .RDS, for now let's keep them there

##Let's place functions that do not rely on the server here, e.g.
#########################Codon usage related###############################################
##Sample a codon
sampcod=function(aa,list,cai){
  newcod=sample((list[[aa]])[,6],1,prob=(list[[aa]])[,cai])
  return(toupper(as.character(newcod)))
}

##Construct new sequence based on multiple samplings
repcds=function(x,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  nnseq=c()
  for(i in seq(1,length(vecseq),by=3)){
    nncod=sampcod(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),list,cai)
    nnseq=append(nnseq,nncod)
  }
  return(paste(nnseq,sep="",collapse=""))
}

##Sample a codon different from input
sampnewcod=function(aa,oldcodon,list,cai){
  oldcodon=toupper(oldcodon)
  if(nrow(list[[aa]]) < 2 ){return(oldcodon)}
  oldcodon =which(as.character(rownames(list[[aa]])) == oldcodon)
  newcod=sample((list[[aa]])[-c(oldcodon),6],1,prob=(list[[aa]])[-c(oldcodon),cai])
  return(toupper(as.character(newcod)))
}

##Construct new sequences without repeating old codons
repnewcds=function(x,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  nnseq=c()
  for(i in seq(1,length(vecseq),by=3)){
    nncod=sampnewcod(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[i:(i+2)]),sep="",collapse=""),list,cai)
    nnseq=append(nnseq,nncod)
  }
  return(paste(nnseq,sep="",collapse=""))
}

##Modify particular positions
modbyposiz=function(x,starts,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  if(length(starts)>0){starts=unique(starts - (starts %% 3) +1 )}
  for(pos in c(starts)){
    nncod=sampnewcod(as.character(tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),list,cai)
    vecseq[pos:(pos+2)]=unlist(strsplit(nncod,""))
  }
  return(paste(vecseq,sep="",collapse=""))
}

##################Related to piRNA search#################################################
##Pis
countpies=function(x,y){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) != 21){return(c())}
  return((stringdist(paste(vecseq[1:14],collapse=""),unlist(y[paste(vecseq[15:20],collapse="")]), method="h")))
}

Strcountpies=function(x,y){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) < 21){return(c())}
  con=0
  for(i in 1:(length(vecseq)-20)){
    con=con+countpies(paste(vecseq[i:(i+20)],collapse=""),y)
  }
  return(con)
}

countmatpies=function(x,y){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) != 21){return(c())}
  return(length(which(stringdist(paste(vecseq[1:14],collapse=""),unlist(y[paste(vecseq[15:20],collapse="")]), method="h")<6)))
}

Strcountmatpies=function(x,y){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) < 21){return(c())}
  con=0
  for(i in 1:(length(vecseq)-20)){
    con=con+countmatpies(paste(vecseq[i:(i+20)],collapse=""),y)
  }
  return(con)
}

condepies=function(x,pies,mm){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) != 20){return(c())}
  return(sum(stringdist(x,pies,method="hamming") <= mm))
}

Strcondepies=function(x,y,mm){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) < 20){return(c())}
  con=0
  for(i in 1:(length(vecseq)-19)){
    con=con+condepies(paste(vecseq[i:(i+19)],collapse=""),y,mm)
  }
  return(con)
}

findpies=function(x,pies,mm){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) != 20){return(c())}
  idx=c(stringdist(x,pies,method="hamming") <= mm)
  if(sum(idx)>0){return(pies[idx])}else{return()}
}

Strfindpies=function(x,y,mm){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) < 20){return(c())}
  con=c()
  for(i in 1:(length(vecseq)-19)){
    con=append(con,findpies(paste(vecseq[i:(i+19)],collapse=""),y,mm))
  }
  return(con)
}

###############################Related to HTML production and sequence viewer###################
###HTML Gene info
CalculateGC = function(x){
  if(!is.character(x)){return(c())}
  x=toupper(x)
  vecseq=unlist(strsplit(x,""))
  return((countPattern("C",x)+countPattern("G",x))/length(vecseq))
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

CalculateCAI = function(x,tabibi,cai,list){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  if(cai == 5 ){cai=11}
  CAIvalues=c()
  for(pos in seq(1,length(vecseq),by=3)){
    am=tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]
    tabs=list[[as.character(am)]]
    CAIvalues=append(CAIvalues,tabs[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),cai]/max(tabs[,cai]))
  }
  return(gm_mean(CAIvalues))
}

HigPrintSeq= function(x,tabpatcol){
  if(!is.character(x)){return(c())}
  if(is.null(nrow(tabpatcol))){return(paste(x))}
  x=toupper(x)
  vecseq=unlist(strsplit(x,""))
  for(i in 1:nrow(tabpatcol)){
    tag=paste0("<span style=\"background-color:",tabpatcol[i,3],"\">")
    vecseq[as.integer(tabpatcol[i,1])]=paste(tag,vecseq[as.integer(tabpatcol[i,1])],sep="")
    vecseq[as.integer(tabpatcol[i,2])]=paste(vecseq[as.integer(tabpatcol[i,2])],"</span>",sep="")
  }
  return(paste(vecseq,sep="",collapse=""))
}

CDSHTMLinfo = function(x,tabibi,cai,list,tabpatcol){
  if (is.null(x)) return(NULL)
  x=toupper(x)
  if(!is.character(x)){return(c())}
  
  if(is.null(nrow(tabpatcol))){colseq=x}else{
    colseq=HigPrintSeq(x,tabpatcol)
    
    ann=c()
    paco=unique(tabpatcol[,c(3,4)])
    for(j in 1:nrow(paco)){
      ann=append(ann,paste0("<span style=\"background-color:",paco[j,1],"\">",paco[j,2],"</span><br>"))
    }
    colseq=paste(c(ann,colseq),sep="",collapse="")
  }
  
  paste0("<b>Codon Adaptation Index</b> based on codon usage selected: ", CalculateCAI(x,tabibi,cai,list), 
         "<br><b>GC content</b>: ", as.integer((CalculateGC(x))*100), "%<br>",
         "<p align=\"justify\"><tt>",
         colseq,
         "</tt></p>"
  )
}

GCHTMLinfo = function(x,tabibi,cai,list){
  if (is.null(x)) return(NULL)
  x=toupper(x)
  if(!is.character(x)){return(c())}
  
  paste0("<b>Bringman Codon Adaptation Index scorex</b>: ", round(CalculateCAI(x,tabibi,cai,list),2), 
         "<br><b>GC content</b>: ", as.integer((CalculateGC(x))*100), "%<br>"
  )
}

SequenceViewer = function(title,div_id,sequence,patterns,colors,tooltips){
  if (is.null(sequence)){return(NULL)}
  if (is.null(div_id)){return(NULL)}
  if(!is.character(sequence)){return(c())}
  if(length(patterns) != length(colors)){return(c())}
  if(length(tooltips) != length(colors)){return(c())}
  
  if(length(patterns) != length(unique(patterns))){return(c())}
  
  if(length(patterns)==0){
    paste0("<div id=\"",div_id,"\"/></div>",
           "<script type=\"text/javascript\">",
           "var seq",div_id," = new Sequence(\'",
           sequence,
           "\');",
           "seq",div_id,".render(\'#",div_id,"\',{",
           "\'showLineNumbers\': true,
  \'wrapAminoAcids\': true,
  \'charsPerLine\': 100,
  \'toolbar\': false,
  \'search\': false,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
</script>                   
                   ")
  }else{
    patitos=matchPattern(DNAString(paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[1:3],sep="",collapse="")),DNAString(paste(sequence,sep="",collapse="")),fixed=T)[1]
    #patitos=1
    for(pat in patterns){
      patitos=c(patitos,matchPattern(DNAString(as.character(pat)),DNAString(paste(sequence,sep="",collapse="")),fixed=T))
    }
    
    
    #subnames=c("Start")
    #subset=patitos[1]
    subnames=c("Start")
    subset=patitos[1]
    
    patitos=patitos[order(start(patitos)),]
    
    if(length(patitos)>1){
      for(n in 2:length(patitos)){
        if(length(disjoin(c(subset,patitos[n]))) != (length(subset)+1)){
          
          if(length(disjoin(c(subset,patitos[n]))) == (length(subset)+2)){
            subset=disjoin(c(subset,patitos[n]))
            subnames=c(subnames,as.character(patitos[n]), paste(as.character(patitos[n-1]),"_;_",as.character(patitos[n]),sep=""))
          }
          
          if(length(disjoin(c(subset,patitos[n]))) == (length(subset))){
            subnames[n-1]=paste(as.character(subnames[n-1]),"_;_",as.character(patitos[n]),sep="")
          }
          
        }else{
          subset=c(subset,patitos[n])
          subnames=c(subnames,as.character(patitos[n]))
          
        }     
        
      }}
    
    #subcol=c("green")
    #subtol=c("ATG")
    
    subcol=c("green")
    subtol=c("Start")
    
    if(length(subnames)>1){
      for(n in 2:length(subnames)){
        if(subnames[n] %in% patterns){
          subcol[n]=colors[which(subnames[n]==patterns)]
          subtol[n]=tooltips[which(subnames[n]==patterns)]
        }else{
          subcol[n]="grey"
          #subtol[n]=subnames[n]
          vecna=unlist(strsplit(subnames[n],"_;_"))
          compna=c()
          for(ja in vecna){
            compna=append(compna,tooltips[which(ja==patterns)])
          }
          subtol[n]=paste(compna,collapse=";")
        }
      }
    }
    
    seqcoverage="{}"
    stpos=start(subset)
    edpos=end(subset)
    if(length(stpos)>0){
      for(s in 1:length(stpos)){
        seqcoverage=paste(seqcoverage,",
                                    {start: ",stpos[s]-1,", end: ",edpos[s],", color: \"white\", bgcolor: \"",subcol[s],"\", underscore: false, tooltip: \"",subtol[s],"\"}",sep="",collapse="")
      }
    }
    
    paste0("<div id=\"",div_id,"\"/></div>",
           "<script type=\"text/javascript\">",
           "var seq",div_id," = new Sequence(\'",
           sequence,
           "\');",
           "seq",div_id,".render(\'#",div_id,"\',{",
           "\'showLineNumbers\': true,
  \'wrapAminoAcids\': true,
  \'charsPerLine\': 100,
  \'toolbar\': false,
  \'search\': false,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
                 ",
"var Sequence",div_id,"Coverage = [",seqcoverage,"];",

"seq",div_id,".coverage(Sequence",div_id,"Coverage);",

"</script>")
    
  }
}

#####Modify function to make it work better
NewSequenceViewer = function(title,div_id,sequence,patterns,colors,tooltips,dflegends,starseq="",endseq=""){
  if (is.null(sequence)){return(NULL)}
  if (is.null(div_id)){return(NULL)}
  if(!is.character(sequence)){return(c())}
  if(length(patterns) != length(colors)){return(c())}
  if(length(tooltips) != length(colors)){return(c())}
  
  if(length(patterns) != length(unique(patterns))){return(c())}
  
  if(starseq==""){starseq=paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[1:10],sep="",collapse="")}
  if(endseq==""){endseq=paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[(nchar(sequence)-10):nchar(sequence)],sep="",collapse="")}
  
  if(length(patterns)==0){
    paste0("<div id=\"",div_id,"\"/></div>",
           "<script type=\"text/javascript\">",
           "var seq",div_id," = new Sequence(\'",
           sequence,
           "\');",
           "seq",div_id,".render(\'#",div_id,"\',{",
           "\'showLineNumbers\': true,
  \'wrapAminoAcids\': true,
  \'charsPerLine\': 100,
  \'toolbar\': false,
  \'search\': false,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
</script>                   
                   ")
  }else{
    ##New paradigm, just add start and endseq to patterns
    patterns=append(patterns,starseq)
    colors=append(colors,"green")
    tooltips=append(tooltips,"Start")
    
    patterns=append(patterns,endseq)
    colors=append(colors,"red")
    tooltips=append(tooltips,"End")
    
    patitos=c()
    
    for(pat in patterns){
      patitos=append(patitos,matchPattern(DNAString(as.character(pat)),DNAString(paste(sequence,sep="",collapse="")),fixed=T))
    }
    
    patitos=patitos[order(start(patitos)),]
    
    subset=patitos[1]
    subnames=as.character(patitos[1])
    
    
    if(length(patitos)>1){
      for(n in 2:length(patitos)){
        if(length(disjoin(c(subset,patitos[n]))) != (length(subset)+1)){
          
          if(length(disjoin(c(subset,patitos[n]))) == (length(subset)+2)){
            subset=disjoin(c(subset,patitos[n]))
            subnames=c(subnames,as.character(patitos[n]), paste(as.character(patitos[n-1]),"_;_",as.character(patitos[n]),sep=""))
          }
          
          if(length(disjoin(c(subset,patitos[n]))) == (length(subset))){
            subnames[n-1]=paste(as.character(subnames[n-1]),"_;_",as.character(patitos[n]),sep="")
          }
          
        }else{
          subset=c(subset,patitos[n])
          subnames=c(subnames,as.character(patitos[n]))
          
        }     
        
      }}
    
    subcol=c("")
    subtol=c("")
    
    if(length(subnames)>1){
      for(n in 1:length(subnames)){
        if(subnames[n] %in% patterns){
          subcol[n]=colors[which(subnames[n]==patterns)]
          subtol[n]=tooltips[which(subnames[n]==patterns)]
        }else{
          subcol[n]="grey"
          #subtol[n]=subnames[n]
          vecna=unlist(strsplit(subnames[n],"_;_"))
          compna=c()
          for(ja in vecna){
            compna=append(compna,tooltips[which(ja==patterns)])
          }
          subtol[n]=paste(compna,collapse=";")
        }
      }
    }
    
    seqcoverage="{}"
    stpos=start(subset)
    edpos=end(subset)
    if(length(stpos)>0){
      for(s in 1:length(stpos)){
        seqcoverage=paste(seqcoverage,",
                                    {start: ",stpos[s]-1,", end: ",edpos[s],", color: \"white\", bgcolor: \"",subcol[s],"\", underscore: false, tooltip: \"",subtol[s],"\"}",sep="",collapse="")
      }
    }
    
    
    newdf=data.frame(Type=c("Start","End","Multiple annotations"),Color=c("green","red","grey"))
    
    
    if(is.data.frame(dflegends)){
      newdf=rbind(newdf,dflegends)
      }
    
    newdf=unique(newdf)
    
    LegendList=paste("{name: \"",newdf[1,1],"\", color: \"",newdf[1,2],"\", underscore: false}",sep="",collapse="")
    for(n in 2:nrow(newdf)){
      LegendList=paste(LegendList,",
                                    {name: \"",newdf[n,1],"\", color: \"",newdf[n,2],"\", underscore: false}",sep="",collapse="")
      }
    
    
    paste0("<div id=\"",div_id,"\"/></div>",
           "<script type=\"text/javascript\">",
           "var seq",div_id," = new Sequence(\'",
           sequence,
           "\');",
           "seq",div_id,".render(\'#",div_id,"\',{",
           "\'showLineNumbers\': true,
  \'wrapAminoAcids\': true,
  \'charsPerLine\': 100,
  \'toolbar\': false,
  \'search\': false,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
                 ",
"var Sequence",div_id,"Coverage = [",seqcoverage,"];

",

"var Sequence",div_id,"Legend = [
",LegendList,"
];

",
"seq",div_id,".coverage(Sequence",div_id,"Coverage);
",

"seq",div_id,".addLegend(Sequence",div_id,"Legend);

",

"</script>")
    
  }
}

############################################Make Ape file#################################
PasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,tabibi,cai,list,PiesList,extracomments=c()){
  if(is.null(sequence)){return(NULL)}
  if(!is.character(sequence)){return(c())}
  if(length(patterns) < 1 ){return(c(paste(sequence)))}
  if(length(patterns) != length(FWDcolors)){return(c())}
  if(length(REVcolors) != length(FWDcolors)){return(c())}
  if(length(tooltips) != length(FWDcolors)){return(c())}
  
  CAIS=CalculateCAI(sequence,tabibi,cai,list) 
  GCp=as.integer((CalculateGC(sequence))*100)
  NoPies=length(Strfindpies(sequence,PiesList,4))
  
  ##Save Lines
  FileLines=c()
  FileLines=append(FileLines,paste("LOCUS",paste(locus_name,sep="",collapse=""),paste(nchar(sequence),"bp ds-DNA", sep=""),"linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "))
  FileLines=append(FileLines,paste("DEFINITION",".",sep="     "))
  FileLines=append(FileLines,paste("ACCESSION",".",sep="     "))
  FileLines=append(FileLines,paste("VERSION",".",sep="     "))
  FileLines=append(FileLines,paste("SOURCE",".",sep="     "))
  FileLines=append(FileLines,paste("ORGANISM","C.elegans",sep="     "))
  
  FileLines=append(FileLines,paste("COMMENT",paste(locus_name),sep="     "))
  
  FileLines=append(FileLines,paste("COMMENT",paste("Codon Adaptation Index",as.character(CAIS)),sep="     "))
  FileLines=append(FileLines,paste("COMMENT",paste("%GC",as.character(GCp)),sep="     "))
  FileLines=append(FileLines,paste("COMMENT",paste("piRNAs sites with less than 4 mismatches seen in this sequence",as.character(NoPies)),sep="     "))
  
  if(length(extracomments)>0){
    for(comocomo in extracomments){
    FileLines=append(FileLines,paste("COMMENT",paste(paste(comocomo,sep="",collapse="")),sep="     "))
    }
    }
  
  posipat=c()
  ##Match sequences
  for(i in 1:length(patterns)){
    stpos=start(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
    edpos=end(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
    if(length(stpos)>0){
      posipat=rbind(posipat, cbind(stpos,edpos,rep(tooltips[i],length(stpos)),rep(FWDcolors[i],length(stpos)),rep(REVcolors[i],length(stpos))))
    }
  }
  
  if(!(is.null(posipat))){
    colnames(posipat)=c("start","end","label","fwdc","revc")
  }
  
  if(!(is.null(posipat))){
    for(i in 1:length(patterns)){
      FileLines=append(FileLines,paste("COMMENT",paste(as.character(tooltips[i]),as.character(patterns[i])),sep="     "))
    }
  }
  
  FileLines=append(FileLines,paste("COMMENT","Generated using wormbuilder.dev",sep="     "))
  FileLines=append(FileLines,paste("COMMENT","ApEinfo:methylated:1",sep="     "))
  
  if(!(is.null(posipat))){
    FileLines=append(FileLines,paste("FEATURES             Location/Qualifiers",sep=""))
    for(n in 1:nrow(posipat)){
      FileLines=append(FileLines,paste("     primer_bind     ",c(posipat[n,1]),"..",c(posipat[n,2]),"",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"",c(posipat[n,4]),"\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"",c(posipat[n,5]),"\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
    }
    
  }
  
  FileLines=append(FileLines,paste("ORIGIN"))
  
  Compseq=unlist(strsplit(sequence,""))
  
  partseq=c()
  
  for(i in seq(1,length(Compseq),10)){
    endseq=i+9
    if(length(Compseq)-i < 9){endseq=length(Compseq)}
    partseq=append(partseq,paste(Compseq[i:endseq],collapse=""))
    
  }
  
  i=1
  for(num in seq(1,length(Compseq),60)){
    index=as.character(num)
    spaces=paste(rep(" ",6-nchar(index)),collapse="")
    endseq=i+5
    if((length(partseq)-i) < 5){endseq=length(partseq)}
    FileLines=append(FileLines , paste(spaces,index," ",paste(partseq[i:(endseq)],collapse=" "),sep=""))
    
    i=i+6
  }
  
  FileLines=append(FileLines,paste("//"))
  
  return(FileLines)
}

##############Twist checking
LogicRange=function(x,treshold){
  if(!is.logical(x)){return(c())}
  flag=FALSE
  st=c()
  ed=c()
  for(idx in 1:length(x)){
    if(flag){
      if(!x[idx]){
        flag=FALSE
        ed=append(ed,idx-1)
      }
    }else{
      if(x[idx]){
        st=append(st,idx)
        flag=TRUE
      }
    }
  }
  tab=c()
  if(length(st) != length(ed)){ed=append(ed,idx)}
  if(length(st) > 0){
    tab=cbind(start=st,end=ed,dist=((ed-st)+1))
  }
  res=tab[which(tab[,"dist"] >= treshold),]
  return(res)
}

LogicHomoRange=function(x,treshold){
  if(!is.logical(x)){return(c())}
  flag=FALSE
  st=c()
  ed=c()
  for(idx in 1:length(x)){
    if(flag){
      if(!x[idx]){
        flag=FALSE
        ed=append(ed,idx)
      }
    }else{
      if(x[idx]){
        st=append(st,idx)
        flag=TRUE
      }
    }
  }
  tab=c()
  if(length(st) != length(ed)){ed=append(ed,idx)}
  if(length(st) > 0){
    tab=cbind(start=st,end=ed,dist=((ed-st)+1))
  }
  res=tab[which(tab[,"dist"] >= treshold),]
  return(res)
}

MicroHomPeWin=function(seq,window){
  if(!is.character(seq)){return(c())}
  if(!is.numeric(window)){return(c())}
  
  seqs=SplitSePeWin(seq,window)
  homos=c()
  for(i in 1:(length(seqs)-window)){
    if(seqs[i] == seqs[i+window]){
      homos=append(homos,TRUE)
    }else{
      homos=append(homos,FALSE)  
    }
  }
  return(homos)
}

SplitSePeWin = function(seq,window){
  if(!is.character(seq)){return(c())}
  if(length(seq) > 1){return(c())}
  if(nchar(seq) < window){return(c(seq))}
  if(window < 1){return(c(seq))}
  
  seq=toupper(seq)
  vecseq=unlist(strsplit(seq,""))
  
  seqs=c()
  for(i in 1:(length(vecseq) - window + 1)){
    seqs=append(seqs, paste(vecseq[i:(i+window-1)],sep="",collapse=""))
  }
  return(seqs)
}

DupSeqsWin=function(seq,window){
  seqs=SplitSePeWin(seq,window)
  revseqs=paste(reverseComplement(DNAStringSet(seqs)))
  dupidx=duplicated(c(seqs,revseqs))
  dupidx=dupidx[1:length(seqs)]
  res=c()
  if(sum(dupidx) > 0){res=unique(seqs[c(dupidx)])}
  return(res)
}

CalculateTM = function(x){
  if(!is.character(x)){return(c())}
  x=toupper(x)
  vecseq=unlist(strsplit(x,""))
  return(4*(countPattern("C",x)+countPattern("G",x)) + 2*(countPattern("A",x)+countPattern("T",x)))
}

TMWindow = function(seq, window){
  if(!is.character(seq)){return(c())}
  if(length(seq) > 1){return(c())}
  if(nchar(seq) < window){return(c())}
  if(window < 1){return(c())}
  
  seq=toupper(seq)
  vecseq=unlist(strsplit(seq,""))
  
  sts=c()
  gcs=c()
  for(i in 1:(length(vecseq) - window + 1)){
    tgc=CalculateTM(paste(vecseq[i:(i+window-1)],sep="",collapse=""))
    sts=append(sts, i)
    gcs=append(gcs, tgc)
  }
  dt=cbind(start=sts,end=(sts+window-1), tm=gcs)
  return(dt)
}

GCWindow = function(seq, window){
  if(!is.character(seq)){return(c())}
  if(length(seq) > 1){return(c())}
  if(nchar(seq) < window){return(c())}
  if(window < 1){return(c())}
  
  seq=toupper(seq)
  vecseq=unlist(strsplit(seq,""))
  
  sts=c()
  gcs=c()
  for(i in 1:(length(vecseq) - window + 1)){
    tgc=CalculateGC(paste(vecseq[i:(i+window-1)],sep="",collapse=""))
    sts=append(sts, i)
    gcs=append(gcs, tgc)
  }
  dt=cbind(start=sts,end=(sts+window-1), gc=gcs)
  return(dt)
}

FlagGCWindow = function(seq, window, loB, hiB){
  if(!is.character(seq)){return(c())}
  if(length(seq) > 1){return(c())}
  if(nchar(seq) < window){return(c())}
  if(window < 1){return(c())}
  seq=toupper(seq)
  vecseq=unlist(strsplit(seq,""))
  
  sts=c()
  eds=c()
  gcs=c()
  
  for(i in 1:(length(vecseq) - window + 1)){
    tgc=CalculateGC(paste(vecseq[i:(i+window-1)],sep="",collapse=""))
    if((tgc > hiB) | (tgc < loB)){
      sts=append(sts, i)
      eds=append(eds, (i+window-1))
      gcs=append(gcs, tgc)
    }
  }
  dt=c()
  if(length(sts)>1){dt=cbind(start=sts, end=eds, gc=gcs)}
  return(dt)
}

##Final Twist function to rule them all
CheckTwistSynthesis = function(sequence){
  message="pass"
  
  #Check if at least is a character
  if(!(is.character(sequence))){
    message = "Not a sequence" 
    return (message)}
  
  ##Now split into characters for DNA sequence
  seqDNA=unlist(strsplit(toupper(sequence),""))
  
  ##Check if DNA sequence
  if((nchar(gsub("A|T|C|G","",toupper(paste(seqDNA,sep="",collapse="")))) != 0)){
    message = "Strange characters in DNA sequence"
    return(message)
  }
  
  ##Check for errors in size
  #lower than 300 bp
  if(length(seqDNA) < 300){ ##Check for errors in size
    message = "Sequence is smaller than 300 bp"
    return(message)
  }
  
  #larger than 500 bp
  if(length(seqDNA) > 5000){ ##Check for errors in size
    message = "Sequence is larger than 5 kb"
    return(message)
  }
  
  ###Now do base composition analysis
  ##GC content
  GCc= CalculateGC(paste(seqDNA,sep="", collapse=""))
  
  if(GCc > .65){
    message = "GC content higher than 65%"
    return(message)
  }
  
  if(GCc < .25){
    message = "GC content lower than 25%"
    return(message)
  }
  
  ##Now compositional analysis per window
  #Check GC content per 50bp window
  GCt= GCWindow(paste(seqDNA,sep="", collapse=""), 50)
  difgc=max(GCt[,"gc"]) - min(GCt[,"gc"])
  
  if(difgc > .52){
    message = "GC content difference between the 50bp windows is larger than 52%"
    return(message)
  }
  
  #Strange double evening out event but seems to be working... what?
  ##I think part of the code of twist is to smooth the gc content on some regions and then re-run analysis
  ##I think, what Amhed in the past did was to smooth by 50 windows and the check once again complex regions
  
  avrm=c()
  for(i in 1:(nrow(GCt)-49)){avrm=append(avrm,mean(GCt[i:(i+49),3]))}
  
  ##Smoth region of GC content
  if(min(avrm) < .2){
    message = "Complex region with low levels of GC detected"
    return(message)
  }
  
  if(max(avrm) > .8){
    message = "Complex region with high levels of GC detected"
    return(message)
  }
  
  ##Now homopolymer track analysis
  ##20bp repeated sequences
  ##First we identify duplicated sequences in the full DNA sequence, by chopping into 20 mers and see if these are present once again in the sequence
  dupseqs=DupSeqsWin(paste(seqDNA,sep="", collapse=""),20)
  
  #COunt how many duplicated sequences
  if(length(dupseqs) > 0){
    #Aerr=append(Aerr,paste("Duplicated 20-mer:",dupseqs))
    #ErrorFlag= ErrorFlag +1
    message = "At least there are 2 20-mer sequences duplicated in the sequence"
    return(message)
  }
  
  ##Now melting temperature for those regions
  ##Calculate 20-mers with TM higher than 60
  TMt=TMWindow(paste(seqDNA,sep="", collapse=""),20)
  timd=which(TMt[,"tm"] > 80)
  
  #Check if the exist
  if(length(timd)>0){
    #Aerr=append(Aerr,paste("20bp regions with high TM (> 80C):",length(timd)))
    #ErrorFlag= ErrorFlag +1
    message = paste("20bp regions with high TM (> 80C):",length(timd))
    return(message)
  }
  
  ##Highest microhomologies on kmer analysis
  kmers=c()
  for(k in 1:20){
    kmers=append(kmers,sum(MicroHomPeWin(paste(seqDNA,sep="", collapse=""),k)))
  }
  
  ##Check
  if(which.max(kmers) > 10){
    #Aerr=append(Aerr,paste("Most frequent micro-homologies seen at ",which.max(kmers), "bps",sep="", collapse=""))
    #ErrorFlag= ErrorFlag +1
    message = paste("Most frequent micro-homologies seen at ",which.max(kmers), "bps",sep="", collapse="")
    return(message)
  }
  
  ##Homo-polymer track larger than 10bp
  t1s=c(MicroHomPeWin(paste(seqDNA,sep="", collapse=""),1),FALSE)
  homotra=LogicHomoRange(t1s,10)
  
  ##Homopolymer tracks, finally
  if(length(homotra)>0){
    #Aerr=append(Aerr,paste("Homo-polymer tracks:",paste(homotra,sep="-", collapse="-")))
    #ErrorFlag= ErrorFlag +1
    message = paste("Homo-polymer tracks:",paste(homotra,sep="-", collapse="-"))
    return(message)
  }
  
  return(message)	
}

##Main function
TwistSynthesisWithCoordinates = function(sequence){
  message="pass"
  errors=c()
  regions=data.frame(Start=integer(),End=integer(),Description=character(), Color=character(),stringsAsFactors=FALSE)
  
  #Check if at least is a character
  if(!(is.character(sequence))){
    message = "Not a sequence" 
    errors=append(errors,message)
    #return (message)
  }
  
  ##Now split into characters for DNA sequence
  seqDNA=unlist(strsplit(toupper(sequence),""))
  
  ##Check if DNA sequence
  if((nchar(gsub("A|T|C|G","",toupper(paste(seqDNA,sep="",collapse="")))) != 0)){
    message = "Strange characters in DNA sequence"
    errors=append(errors,message)
    #return (message)
  }
  
  ##Check for errors in size
  #lower than 300 bp
  if(length(seqDNA) < 300){ ##Check for errors in size
    message = "Sequence is smaller than 300 bp"
    errors=append(errors,message)
    #return (message)
  }
  
  #larger than 500 bp
  if(length(seqDNA) > 5000){ ##Check for errors in size
    message = "Sequence is larger than 5 kb"
    errors=append(errors,message)
    #return (message)
  }
  
  ###Now do base composition analysis
  ##GC content
  GCc= CalculateGC(paste(seqDNA,sep="", collapse=""))
  
  if(GCc > .65){
    message = "GC content higher than 65%"
    errors=append(errors,message)
    #return (message)
  }
  
  if(GCc < .25){
    message = "GC content lower than 25%"
    errors=append(errors,message)
    #return (message)
  }
  
  ##Now compositional analysis per window
  #Check GC content per 50bp window
  GCt= GCWindow(paste(seqDNA,sep="", collapse=""), 50)
  difgc=max(GCt[,"gc"]) - min(GCt[,"gc"])
  
  if(difgc > .52){
    message = "GC content difference between the 50bp windows is larger than 52%"
    regions[nrow(regions)+1,]=c(GCt[which.max(GCt[,"gc"])[1],1],GCt[which.max(GCt[,"gc"])[1],2],"Highest GC 50bp window","purple")
    regions[nrow(regions)+1,]=c(GCt[which.min(GCt[,"gc"])[1],1],GCt[which.min(GCt[,"gc"])[1],2],"Lowest GC 50bp window","blue")
    errors=append(errors,message)
    #return (message)
  }
  
  #Strange double evening out event but seems to be working... what?
  ##I think part of the code of twist is to smooth the gc content on some regions and then re-run analysis
  ##I think, what Amhed in the past did was to smooth by 50 windows and the check once again complex regions
  
  avrm=c()
  for(i in 1:(nrow(GCt)-49)){avrm=append(avrm,mean(GCt[i:(i+49),3]))}
  
  ##Smooth region of GC content
  if(min(avrm) < .2){
    regions[nrow(regions)+1,]=c(GCt[which.min(avrm)[1],1],GCt[which.min(avrm)[1]+49,2],"Complex region with low levels of GC","green")
    message = "Complex region with low levels of GC detected"
    errors=append(errors,message)
    #return (message)
  }
  
  if(max(avrm) > .8){
    regions[nrow(regions)+1,]=c(GCt[which.max(avrm)[1],1],GCt[which.max(avrm)[1]+49,2],"Complex region with high levels of GC","pink")
    message = "Complex region with high levels of GC detected"
    errors=append(errors,message)
    #return (message)
  }
  
  ##Now homopolymer track analysis
  ##20bp repeated sequences
  ##First we identify duplicated sequences in the full DNA sequence, by chopping into 20 mers and see if these are present once again in the sequence
  dupseqs=DupSeqsWin(paste(seqDNA,sep="", collapse=""),20)
  
  #COunt how many duplicated sequences
  if(length(dupseqs) > 0){
    #Aerr=append(Aerr,paste("Duplicated 20-mer:",dupseqs))
    #ErrorFlag= ErrorFlag +1
    unidupseqs=unique(dupseqs)
    for(pat in unidupseqs){
      matches=matchPattern(DNAString(pat),DNAString(paste(seqDNA,sep="", collapse="")),fixed=T)
      regions[(nrow(regions)+1):(nrow(regions)+length(matches)),]=c(start(matches),end(matches),paste("Duplicated sequence:",as.character(matches)),rep("gold",length(matches)))
    }
    message = paste("There are at least",length(dupseqs), "20-mer sequences identical in DNA sequence")
    errors=append(errors,message)
    #return (message)
  }
  
  ##Now melting temperature for those regions
  ##Calculate 20-mers with TM higher than 60
  TMt=TMWindow(paste(seqDNA,sep="", collapse=""),20)
  timd=which(TMt[,"tm"] > 80)
  
  #Check if the exist
  if(length(timd)>0){
    #Aerr=append(Aerr,paste("20bp regions with high TM (> 80C):",length(timd)))
    #ErrorFlag= ErrorFlag +1
    regions[(nrow(regions)+1):(nrow(regions)+length(timd)),]=c(TMt[timd,1],TMt[timd,2],rep("Region with high melting temperature",length(timd)),rep("purple",length(timd)))
    message = paste("20bp regions with high TM (> 80C):",length(timd))
    errors=append(errors,message)
    #return (message)
  }
  
  ##Highest microhomologies on kmer analysis
  kmers=c()
  for(k in 1:20){
    kmers=append(kmers,sum(MicroHomPeWin(paste(seqDNA,sep="", collapse=""),k)))
  }
  
  ##Check
  if(which.max(kmers) > 10){
    #Aerr=append(Aerr,paste("Most frequent micro-homologies seen at ",which.max(kmers), "bps",sep="", collapse=""))
    #ErrorFlag= ErrorFlag +1
    message = paste("Most frequent micro-homologies seen at ",which.max(kmers), "bps",sep="", collapse="")
    errors=append(errors,message)
    #return (message)
  }
  
  ##Homo-polymer track larger than 10bp
  t1s=c(MicroHomPeWin(paste(seqDNA,sep="", collapse=""),1),FALSE)
  homotra=LogicHomoRange(t1s,10)
  
  ##Homopolymer tracks, finally
  if(length(homotra)>0){
    #Aerr=append(Aerr,paste("Homo-polymer tracks:",paste(homotra,sep="-", collapse="-")))
    #ErrorFlag= ErrorFlag +1
    regions[(nrow(regions)+1):(nrow(regions)+nrow(homotra)),]=c(homotra[,1],homotra[,2],rep("Micro-homology found between adjacent regions",nrow(homotra)),rep("cyan",nrow(homotra)))
    message = paste("Homo-polymer tracks:",paste(homotra,sep="-", collapse="-"))
    errors=append(errors,message)
    #return (message)
  }
  
  return(list(message,errors,regions))	
}

##############################################################################SERVER FUNCTION######################################
shinyServer(function(input, output, session) {
  ###########################################################################################
  #########################Related to user session###########################################
  ###########################################################################################
  
  ##Retrieve unique ID for the session
  session_id <- session$token
  ##Create temporary folder for unique user
  system(paste("mkdir -p DATA/users/",session_id,sep=""))
  ###On exit, force the remove of directory
  ##Attention: Interactive sessions create problems because the listening of the server stops in main directory and not in sub directory
  session$onSessionEnded(function(){
    system(paste("rm -rf DATA/users/",session_id,sep=""))
  }
  )
  
  ###Create path so life becomes easier
  UserPath = paste("DATA/users/",session_id,"/",sep="")
  
  ###########################################################################################
  #########################Functions#########################################################
  ###########################################################################################
  ##Start by simple loading
  #output$DynamicUserInterface <- renderUI({HTML("<b>Loading interface</b>")})
  #output$DynamicUserInterface <- renderUI({uiOutput("AllResults")})
  #output$AllResults <- renderUI({fluidRow(uiOutput("button4later"))})
    
  #output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
  
  
  ###########################################################################################
  ##############################Download handlers############################################
  ###########################################################################################
  
  ##################Download handler
  # output$DownSeqOut <- downloadHandler(
  #   filename <- function() {
  #     paste("Genebuild", "fasta", sep=".")
  #   },
  #   
  #   content <- function(file) {
  #     file.copy(paste("DATA/users/",session_id,"/SeqOpop.fasta", sep=""), file)
  #   },
  # )
  # 
  # output$DownBlockConstruct <- downloadHandler(
  #   filename <- function() {
  #     paste("Genebuild", "fasta", sep=".")
  #   },
  #   
  #   content <- function(file) {
  #     file.copy(paste("DATA/users/",session_id,"/BlockConstruct.fasta", sep=""), file)
  #   },
  # )
  # 
  # 
  # output$DownOriApe <- downloadHandler(
  #   filename <- function() {
  #     paste("Input_sequence", "gb", sep=".")
  #   },
  #   
  #   content <- function(file) {
  #     file.copy(paste("DATA/users/",session_id,"/Seqog.gb", sep=""), file)
  #   },
  # )
  # 
  # output$DownOptiApe <- downloadHandler(
  #   filename <- function() {
  #     paste("Output_sequence", "gb", sep=".")
  #   },
  #   
  #   content <- function(file) {
  #     file.copy(paste("DATA/users/",session_id,"/Seqpop.gb", sep=""), file)
  #   },
  # )
  
  
  ###Testing with input name
  output$DownSeqOut <- downloadHandler(
    filename <- function() {
      SeqNameIn=as.character(input$nameinput)
      
      ###Workout name
      SeqNameIn=gsub(" ","", SeqNameIn)
      if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
      
      paste(SeqNameIn, "fasta", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/SeqOpop.fasta", sep=""), file)
    }
  )
  
  output$DownBlockConstruct <- downloadHandler(
    filename <- function() {
      SeqNameIn=as.character(input$nameinput)
      
      ###Workout name
      SeqNameIn=gsub(" ","", SeqNameIn)
      if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
      
      paste(SeqNameIn, "fasta", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/BlockConstruct.fasta", sep=""), file)
    }
  )
  
  
  output$DownOriApe <- downloadHandler(
    filename <- function() {
      SeqNameIn=as.character(input$nameinput)
      
      ###Workout name
      SeqNameIn=gsub(" ","", SeqNameIn)
      if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
      
      paste(SeqNameIn, "gb", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/Seqog.gb", sep=""), file)
    }
  )
  
  output$DownOptiApe <- downloadHandler(
    filename <- function() {
      SeqNameIn=as.character(input$nameinput)
      
      ###Workout name
      SeqNameIn=gsub(" ","", SeqNameIn)
      if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
      
      paste(paste("Optimized",SeqNameIn,sep="_"), "gb", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/Seqpop.gb", sep=""), file)
    }
  )
  

  ###########################################################################################
  #########################Events############################################################
  ###########################################################################################
  
  ###Control panels##########################################################################################
  observeEvent(input$link_to_tabpanel_sequenceadaptation, {
    newvalue <- "Sequence Adaptation"
    updateTabsetPanel(session, "panels", newvalue)
  })
  
  #### UI control ###########
  observeEvent(input$checkMetSites,{
    genz=as.character(input$Genzymes)
    oenz=as.character(input$Oenzymes)
    
    if(input$checkMetSites){
      updatePrettyCheckboxGroup(
        session = session,
        inputId = "Genzymes",
        choices = c(GoldenE[-c(which(GoldenE$DamDcm == "TRUE")),"Enzyme"]), inline=TRUE, selected=genz)
      updatePrettyCheckboxGroup(
        session = session,
        inputId = "Oenzymes",
        choices = c(OtherE[-c(which(OtherE$DamDcm == "TRUE")),"Enzyme"]), inline=TRUE, selected=oenz)
    }else{
      updatePrettyCheckboxGroup(
        session = session,
        inputId = "Genzymes",
        choices = c(GoldenE[,"Enzyme"]), inline=TRUE, selected=genz)
      updatePrettyCheckboxGroup(
        session = session,
        inputId = "Oenzymes",
        choices = c(OtherE[,"Enzyme"]), inline=TRUE,selected=oenz)
    }
    
  })
  
  ##########################################################################################################
  ##############################Main Event that brings the dynamic UI at its initial state##################
  ##########################################################################################################
  ###Render UI controls given a function called init
  inituiui=function(){
    output$DynamicUserInterface <- renderUI({
      fluidRow(
        HTML("<h2><i>C. elegans</i> transgene adaptation</h2>"),
        radioButtons("intypeinput", label = HTML("<h4>Input</h4>"),
                     choices = list("DNA" = 1, 
                                    "Protein" = 2), 
                     selected = 1, inline = TRUE, width='100%'),
        textAreaInput("nameinput", label = HTML(""), value = "", resize="none", placeholder= "Sequence name (optional)", rows=1),
        conditionalPanel(condition = "input.intypeinput==1",
                         fluidRow(
                           column(8,
                                  textAreaInput("seqDNA", label = HTML(""), value = "", cols= 100, rows=5, width = "600px"),
                                  HTML("<h5>Include start and stop codons. Maximum gene length is 10 kb.</h5>")
                           )
                         )),
        conditionalPanel(condition = "input.intypeinput==2",
                         fluidRow(
                           column(8,
                                  textAreaInput("seqPROT", label = HTML(""), value = "", cols= 100, rows=5, width = "600px"),
                                  HTML("Maximum protein length is 3332 aa.<br><br>")
                           ),
                         )),
        
        ###Parameters
        HTML("<h4>Sequence manipulation:</h4>"),
        selectInput("selectCAI", label = HTML("<b>Optimization method
                                                           [<a href=\"\" onclick=\"$('#explain_codon').toggle(); return false;\">info</a>]
                                                           </b>"), 
                    choices = list("--- Codon table ---" = 1, 
                                   "Ubiquitous (stochastic, frequently used codons)" = 2, 
                                   #"Germline" = 3, 
                                   #"Neuronal" = 4, 
                                   #"Somatic" = 5, 
                                   "--- Published algorithms ---" = 6, 
                                   "Redemann et al. Nature Methods (2011)" = 7, 
                                   "Fielmich et al. eLife (2018)" = 8), 
                    selected = 1),
        HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_codon\">
                        Codon frequencies were calculated from the 500 highest expressed tissue-specific genes identified by RNA-seq on cell sorted nuclei. <a href=\"https://doi.org/10.1101/2020.02.20.958579\">Serizay <i>et al.</i> (2020)</a>. <a href=\"https://www.ahringerlab.com\">Ahringer lab</a>.
                          <br> Please note that codon sampling is randomized and the optimized output sequence is therefore not invariant.
                        <br>
                        Max. expression corresponds to Codon Adaptation Index = 1. <a href=\"http://www.nature.com/nmeth/journal/v8/n3/full/nmeth.1565.html\">Redemann <i>et al.</i> (2011)</a>. <a href=\"https://worm.mpi-cbg.de/codons/cgi-bin/optimize.py\"><i>C. elegans</i> Codon Adapter</a>.
                        <br>
                        <u>Ge</u>rm <u>li</u>ne <u>op</u>timization was developed by <a href=\"https://www.utdickinsonlab.org\">Dan Dickinson</a> and was described in <a href=\"https://elifesciences.org/articles/38198\">Fielmich <i>et al.</i> (2018)</a>. See also: <a href=\"http://104.131.81.59/\">http://104.131.81.59/</a>. 
                    </div></p>"),
        
        checkboxInput("checkboxRibo", label = HTML("<b>Optimize ribosomal binding 
                                                                [<a href=\"\" onclick=\"$('#explain_ribo').toggle(); return false;\">info</a>]</b>"), value = FALSE),
        HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_ribo\">
                        This option minimizes the folding energy of positions -4 to +39 to increase ribosomal binding. This option also adds the consensus start sequence (aaaaATG). Note that ribosomal binding site optimization is not compatible with the germline optimization algorithm.
                    </div></p>"),
 fluidRow(
   column(8,
          checkboxInput("checkPirna", label = HTML("<b>Minimize <i>C. elegans</i> piRNA homology
                                                              [<a href=\"\" onclick=\"$('#explain_piRNA').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%')
   )),
 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_piRNA\">
                        This option annotates and minimizes sequence homology to piRNAs to reduce germline silencing. The algorithm removes, when possible, sequences with less than four mismatches to all endogenous piRNAs. Please note that this algorithm is computationally demanding.
                    </div></p>"),
 checkboxInput("checkEnzySites", label = HTML("<b>Remove restriction enzyme sites:
                                                               [<a href=\"\" onclick=\"$('#explain_enzysites').toggle(); return false;\">info</a>]
                                                                   </b>"), value = FALSE, width='100%'),
 conditionalPanel(condition = "input.checkEnzySites==1",
                  fluidRow(
                    #HTML("<tt>"),
                    column(6,
                           prettyCheckboxGroup("Genzymes", "Used for GoldenGate assembly:",
                                               c(GoldenE$Enzyme), inline=TRUE),
                           prettyCheckboxGroup("Oenzymes", "Other:",
                                               c(OtherE$Enzyme), inline=TRUE)
                          # checkboxInput("checkMetSites", label = HTML("<b>Dam/Dcm
                          #                                     [<a href=\"\" onclick=\"$('#explain_metsites').toggle(); return false;\">info</a>]
                         #                                      </b>"), value = FALSE, width='100%'),
                        #   HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_metsites\">
                       # Disable restriction sites affected by <a href=\"https://blog.addgene.org/plasmids-101-methylation-and-restriction-enzymes\">Dam/Dcm</a> methylases. </div></p>")
                    ),
                    #HTML("</tt>")
                  )),
 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_enzysites\">
                        Note, this optimization is performed after piRNA removal and may re-introduce piRNAs sites. Any piRNAs sites are annotated in the sequence output.
                          </div></p>"),
 checkboxInput("checkfouras", label = HTML("<b>Start sequence with consensus start site
                                                              [<a href=\"\" onclick=\"$('#explain_consensusstart').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%'),
 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_consensusstart\">
                        This option also adds the consensus start sequence (aaaaATG).
                    </div></p>"),
 fluidRow(
   column(8,
          checkboxInput("checkIntron", label = HTML("<b>Add three introns
                                                               [<a href=\"\" onclick=\"$('#explain_introns').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
          conditionalPanel(condition = "input.checkIntron==1",
                           radioButtons("intropt", label = HTML(""),
                                        choices = list(
                                          #"Synthetic, Golden Gate compatible (BsaI, 51 bp, 33% GC)" = 1, 
                                          "rps-0 (55 bp, 15% GC)" = 2, 
                                          "rps-5 (65 bp, 22% GC)" = 3,
                                          "rps-20 (62 bp, 28% GC)" = 4,
                                          "Canonical Fire lab introns" = 5
                                        ), 
                                        selected = 2, width='100%'),
                           radioButtons("intdistop",label = HTML("Intron placement"),
                                        choices = list("Early start" = 1, 
                                                       "equidistant" = 2), 
                                        selected = 1, width='100%', inline = TRUE),
                           checkboxInput("checkintframe", label = HTML("Force introns in reading frame"), value = FALSE)
          ))),
 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_introns\">
                        Introns can improve transgene expression (<a href=\"https://pubmed.ncbi.nlm.nih.gov/8244003/\">Okkema <i>et al.,</i> 1993</a>). The introns are indicated by lower-case letters and are inserted at consensus splice sites (<a href=\"https://www.ncbi.nlm.nih.gov/books/NBK20075/\"><i>C. elegans</i> II, 2. ed</a>). 
                        Introns within the first 150 basepairs are particularly efficient at improving germline expression (<a href=\"https://www.nature.com/articles/s41467-020-19898-0\">Al Johani <i>et al.,</i> 2020</a>).
                    </div></p>"),
 ###Addition of UTRs
 fluidRow(
   column(8,
          checkboxInput("checkUTRs", label = HTML("<b>Append UTR's
                                                               [<a href=\"\" onclick=\"$('#explain_UTRs').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
          conditionalPanel(condition = "input.checkUTRs==1",
                           radioButtons("p5UTR", label = HTML("5' UTR"),
                                        choices = list(
                                          "None" = 1, 
                                          "Fire lab synthetic spliced" = 2
                                        ), 
                                        selected = 1, width='100%', inline = TRUE),
                           radioButtons("p3UTR", label = HTML("3' UTR"),
                                        choices = list(
                                          "None" = 1, 
                                          "rps-1 3' UTR" = 2,
                                          "rps-4 3' UTR" = 3,
                                          "tbb-2 3' UTR" = 4
                                        ), 
                                        selected = 1, width='100%', inline = TRUE)
          ))),
 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_UTRs\">
                        Append untranslated regions (UTRs) before and after coding sequences.
                    </div></p>"),
 
 ###Output manipulation
 #HTML("<h4>Visual output:</h4>"),
 #checkboxInput("checkAnno", label = HTML("<b>Annotate sequences
  #                                                            [<a href=\"\" onclick=\"$('#explain_anno').toggle(); return false;\">info</a>]
  #                                                            </b>"), value = FALSE, width='100%'),
 #HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_anno\">
#                        Display the location of piRNA and restriction sites on both input and output sequences.
#                    </div></p>"),
 conditionalPanel(condition = "input.intypeinput==1",
                  checkboxInput("checkAnal", label = HTML("<b>Analytical mode
                                                              [<a href=\"\" onclick=\"$('#explain_anal').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%'),
                  HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_anal\">
                        This mode analyzes the input sequence but does not perform any optimization.
                    </div></p>"),
                  checkboxInput("checkModOnly", label = HTML("<b>Skip codon optimization routine
                                                              [<a href=\"\" onclick=\"$('#explain_modonly').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%'),
                  HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_modonly\">
                        Skip codon optimization but run the rest of the algorithm, e.g. RBS optimization, piRNA removal and addition of extra sequences.
                    </div></p>"),
                  ),
checkboxInput("checkTwisty", label = HTML("<b>Check if sequence can be synthetized into a gene fragment
                                                              [<a href=\"\" onclick=\"$('#explain_Twisty').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%'),
HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_Twisty\">
                        Do not optimize and check if input sequence follows gene design guidelines of <a href=\"https://www.twistbioscience.com/\">Twist biosciences</a>.
                    </div></p>"),
 actionButton("actionSeq", label = "Optimize sequence"),
 verbatimTextOutput("ErrorMessage"),
 hr()
 #uiOutput("AllResults")
      )
    })
    
  }
  
  ##uiOutput("DynamicUserInterface")
  observeEvent(session$clientData,{inituiui()})
  
  ######Render UI controls once reset button is clicked
  observeEvent(input$actionRESET, {
    clearAllResults()
    inituiui()})
  
  ##Clear all results tab
  clearAllResults = function(){
    
    #uiOutput("downloadoptseq"),
    output$downloadoptseq <- renderUI({})
    
    #htmlOutput("oldsequence"),
    output$oldsequence <- renderUI({}) 
    
    #uiOutput("downloadApeOriseq"),
    output$downloadApeOriseq <- renderUI({})
    
    #tableOutput("OriPiTab"),
    output$OriPiTab <- renderTable({})
    
    #htmlOutput("newsequence"),
    output$newsequence <- renderUI({})
    
    #uiOutput("downloadApeOptiseq"),
    output$downloadApeOptiseq <- renderUI({})
    
    #tableOutput("OptiPiTab")
    output$OptiPiTab <- renderTable({})
    
    #outputbutton4later
    #output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
    output$button4later <- renderUI({HTML("<b></b>")})
    }
  
  #######Javascript events
  observeEvent(input$clockactivity, {
    checkfilestatus(UserPath)
  })
  
  checkfilestatus = function(path){
    StatusMessage=""
    if(file.exists(paste(path,"JobStatus.txt",sep=""))){
      StatusMessage=readLines(paste(path,"JobStatus.txt",sep=""))[1]
      }else{
        StatusMessage="Initializing..."
      }
    session$sendCustomMessage("JobStatusMessage", StatusMessage)
    }

  #########################################Main generator after clicking input#######################
  #####Transgene generation#####
  observeEvent(input$actionSeq, {  
    
    ErrorFlag=0
    
    SeqNameIn=as.character(input$nameinput)
    CodonAl=as.integer(input$selectCAI)
    
    ##Check for strange naming, particularly <, >, comas, ., and ;
      if((length(grep(";|>|<",SeqNameIn)) != 0)){ ##Check for strange characters in name
        output$ErrorMessage <- renderText({
          paste("Error: Special characters >, <, and ; are not allowed in the sequence name")
        })
        return(NULL)
      }
      
    #Ori name
    OriSeqNameIn=as.character(input$nameinput)
    if(OriSeqNameIn ==""){OriSeqNameIn="Input sequence"}
    
    ###Workout name
    SeqNameIn=gsub(" ","", SeqNameIn)
    if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
    
    ##Remove piRNas?
    FlaPi=input$checkPirna
    
    ##Add introns?
    FlaIn=input$checkIntron
    
    ##Perform RBS optimization?
    FlaRi=input$checkboxRibo
    
    ##Remove RE sites?
    FlaEnz=input$checkEnzySites
    
    ##Perform only analysis?
    FlaAna=input$checkAnal
    
    ##Trailing Aas?
    Fla5paaa=input$checkfouras
    
    ##Do not codon optimize
    FlaModOnly=input$checkModOnly
    
    ###Patch in wrong way to mantain analysis always active
    #FlaAno=input$checkAnno
    ##Annotate sequences?
    FlaAno=TRUE
    
    typinput=input$intypeinput
    secprousr=gsub(" |\n","", input$seqPROT)
    secdnausr=gsub(" |\n","", input$seqDNA)
    gegegenzymes=input$Genzymes
    gogogonzymes=input$Oenzymes
    
    ##Additional flag so inputs are obtain from the beginning of function
    if(FlaIn){
      typeIn=as.integer(input$intropt)
      typedistint=as.integer(input$intdistop)
      Flaframeint=input$checkintframe
    }
    
    ##Now check introns
    FlaTURS=input$checkUTRs
    befUTR=""
    aftUTR=""
    
    if(FlaTURS){
      befUTR=c("","CCGGGATTGGCCAAAGGACCCAAAGgtatgtttcgaatgatactaacataacatagaacattttcagGAGGACCCTTGGAGGGTACCGGTAG")[as.integer(input$p5UTR)]
      aftUTR=c("","tcttataatttcattgttatgtcgcattgcgataaatgttaaaattaaaaaacttc","actgttgtttttgttgaaaaataaaattgttaatctaaaaa","atgcaaaatcctttcaagcattcccttcttctctatcactcttctttctttttgtcaaaaaattctctcgctaatttatttgcttttttaatgttattattttatgactttttatagtcactgaaaagtttgcatctgagtgaagtgaatgctatcaaaatgtgattct")[as.integer(input$p3UTR)]
      }
    
    if(befUTR != ""){Fla5paaa=TRUE}
    ##CHeck if user wants to look for issues with gene synthesis
    FlaTwisty=input$checkTwisty
  
    output$ErrorMessage <- renderText({""})
    
    output$AllResults <- renderUI({})
    
    seqDNA=""
    
    ##Status file
    statfiletowork = paste(UserPath,"JobStatus.txt",sep="")
    #Remove if present
    if(file.exists(statfiletowork)){file.remove(statfiletowork)}
    
    
    if(as.integer(typinput) == 2){###Input protein
      if((ErrorFlag == 0) & (nchar(gsub("A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y|\\*","",toupper(secprousr))) != 0)){ ##Check for strange non aminoacid characters
        output$ErrorMessage <- renderText({
          paste("Error: Unrecognized characters in sequence:",secprousr)
        })
        ErrorFlag=1
      }else{
        seqPROTO=gsub("*","X",toupper(secprousr),fixed=TRUE)
        trDNA=paste(unlist(sapply(unlist(strsplit(toupper(seqPROTO),"")),function(x){sampcod(x,AAtoCodF,1)})),sep="",collapse="")
        seqDNA=unlist(strsplit(toupper(trDNA),""))
        ###Make sure to make analysis flag off and Twist off
        FlaAna=FALSE
        FlaTwisty=FALSE
      }
    }else{
      seqDNA=unlist(strsplit(toupper(secdnausr),""))
    }
    
    ##Patch for codon algorithm, it should me modified properly later
    if(FlaModOnly){
      CodonAl = 2
    }
    
    ###FLaTwisty
    ##Change error flag so we skip the other checks
    if(FlaTwisty&(ErrorFlag==0)){
      ErrorFlag=5 #5, becuase... why not
    }
    
    ##Twist supercededs analytical mode
    if((FlaTwisty)&(FlaAna)){
      output$ErrorMessage <- renderText({
        paste("Error: Both analysis cannot be active at the same time. Analytical mode will be deactivated")
      })
      FlaAna=FALSE
      updateCheckboxInput(session, "checkAnal", value = FALSE)
      }
    
    ##Analytical mode supercedes codon table
    if(!(FlaAna) & (ErrorFlag == 0)){
      if( (CodonAl == 1) | ( CodonAl == 6)){ ##Check for error on optimization algorithm
        output$ErrorMessage <- renderText({
          paste("Select a codon table or a published algorithm to use for adaptation.")
        })
        ErrorFlag=1
      }else{
        if(CodonAl > 6){CodonAl = CodonAl - 2}else{
          CodonAl = CodonAl - 1
        }
      }
    }
    
    
    if((ErrorFlag == 0) & ((length(seqDNA) %% 3) != 0)){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence is not multiple of three")
      })
      ErrorFlag=1
    }
    
    if((ErrorFlag == 0) & ((length(seqDNA) < 50) & (FlaRi))){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence is too short to optimize for ribosome binding")
      })
      ErrorFlag=1
    }
    
    if((ErrorFlag == 0) & (length(seqDNA) <= 39)){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence length should be at least of 39 bp (13 aa)")
      })
      ErrorFlag=1
    }
    
    if((ErrorFlag == 0) & (length(seqDNA) > 10000)){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence is larger than 10 kb.")
      })
      ErrorFlag=1
    }
    
    if((ErrorFlag == 0) & ((length(seqDNA) < 450) & (FlaIn))){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence is too short to add ~150bp spaced introns. The option will be deactivated")
      })
      FlaIn=FALSE
      updateCheckboxInput(session, "checkIntron", value = FALSE)
    }
    
    if((ErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(paste(seqDNA,sep="",collapse="")))) != 0)){ ##Check for strange non ATCG characters
      output$ErrorMessage <- renderText({
        paste("Error: Unrecognized characters or multiple newlines found in:",paste(seqDNA,sep="",collapse=""))
      })
      ErrorFlag=1
    }
    
    if((ErrorFlag == 0) & (paste(seqDNA[1:3],sep="",collapse="") != "ATG")){ ##Check for errors in size
      if(as.integer(typinput) == 2){
        output$ErrorMessage <- renderText({
          paste("Error: Sequence does not start with Methionine")
        })
      }else{
        output$ErrorMessage <- renderText({
          paste("Error: Sequence does not start with ATG")
        })
      }
      
      ErrorFlag=1
    } 
    
    if((ErrorFlag == 0) & ((FlaPi)&(CodonAl == 5))){ ##Check if CAI = 1 is required
      output$ErrorMessage <- renderText({
        paste("Error: Sequence cannot remove piRNAs and mantain a Codon Adaptation Index equal to 1. piRNAs will not be removed")
      })
      FlaPi=FALSE
      updateCheckboxInput(session, "checkPirna", value = FALSE)
      return(NULL)
    } 
    
    if((ErrorFlag == 0) & ((FlaRi)&(CodonAl == 6))){ ##Check if GLO is required
      output$ErrorMessage <- renderText({
        paste("Error: Sequence cannot be optimized for ribosome binding and Germline expression at the same time. GLO algoritmh will be run on full sequence")
      })
      FlaRi=FALSE
      updateCheckboxInput(session, "checkboxRibo", value = FALSE)
      return(NULL)
    }
    
    if(ErrorFlag == 0){##Check stop codon
      seqiqi=paste(seqDNA,sep="",collapse="")
      stpos=c()
      stpos=start(matchPattern("*",translate(DNAString(seqiqi))))
      if(length(stpos) == 0){
        output$ErrorMessage <- renderText({
          paste("Error: Sequence does not have stop codon")
        })
        ErrorFlag=1
      }
      if(length(stpos) > 1){
        output$ErrorMessage <- renderText({
          paste("Error: Sequence have multiple stop codons")
        })
        ErrorFlag=1
      }
      if(length(stpos) == 1){
        if( (stpos) != (length(translate(DNAString(seqiqi)))) ){
          output$ErrorMessage <- renderText({
            paste("Error: Sequence has a stop codon but not at its end.")
          })
          ErrorFlag=1
        }
      }
    }
    
    ####Passed parameters
    
    eeexxttpar=c()
    if(as.integer(typinput) == 1){eeexxttpar=append(eeexxttpar,paste("Input sequence: ",secdnausr,sep="",collapse=""))}
    if(as.integer(typinput) == 2){eeexxttpar=append(eeexxttpar,paste("Input sequence: ",secprousr,sep="",collapse=""))}
    if(FlaModOnly | FlaAna){eeexxttpar=append(eeexxttpar,"Skip optimization routine: Yes")}else{
      if(CodonAl == 2){eeexxttpar=append(eeexxttpar,"Optimization routine: Ubiquitous (stochastic, frequently used codons)")}
      if(CodonAl == 5){eeexxttpar=append(eeexxttpar,"Optimization routine: Redemann et al. (2011)")}
      if(CodonAl == 6){eeexxttpar=append(eeexxttpar,"Optimization routine: Fielmich et al. (2018)")}
    }
    
    
    if(FlaRi){eeexxttpar=append(eeexxttpar,"RBS optimization: Yes")}else{eeexxttpar=append(eeexxttpar,"RBS optimization: No")}
    if(FlaPi){eeexxttpar=append(eeexxttpar,"Remove piRNA sites: Yes")}else{eeexxttpar=append(eeexxttpar,"Remove piRNA sites: No")}
    if(FlaEnz){eeexxttpar=append(eeexxttpar,"Remove RE sites: Yes")}else{eeexxttpar=append(eeexxttpar,"Remove RE sites: No")}
    if(FlaIn){eeexxttpar=append(eeexxttpar,"Add introns: Yes")}else{eeexxttpar=append(eeexxttpar,"Add introns: No")}
    if(FlaTURS){eeexxttpar=append(eeexxttpar,"Add UTRs: Yes")}else{eeexxttpar=append(eeexxttpar,"Add UTRs: No")}
    if(Fla5paaa){eeexxttpar=append(eeexxttpar,"Add consensus start site: Yes")}else{eeexxttpar=append(eeexxttpar,"consensus start site: No")}
    
    ##Analytical mode supercedes main routine
    ##But in case that analytical mode and piRNAi is off the former routine will be shot down just after
    if((ErrorFlag == 0)&(FlaAna)){
      ##Load all results in dynamic ui
      output$DynamicUserInterface <- renderUI({uiOutput("AllResults")})
      
      output$AllResults <- renderUI({
        fluidRow(
          #actionButton("actionRESET", label = "RESET"),
          uiOutput("button4later"),
          htmlOutput("oldsequence"),
          uiOutput("downloadApeOriseq"),
          hr(),
          tableOutput("OriPiTab")
          #dataTableOutput("OriPiTab")
        )
      })
      
      #output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
      output$button4later <- renderUI({HTML("<b></b>")})
      
      if(as.integer(typinput) == 1){ ###If original sequence
        output$oldsequence <-renderUI({
          seqiqi=toupper(seqiqi)
          
          piss=Strfindpies(seqiqi,Pies,4)
          if(length(piss)>0){
            popos=c()
            papas=c()
            for(pipi in piss){
              patotes=as.character(matchPattern(DNAString(pipi),DNAString(seqiqi),max.mismatch=4,fixed=T))
              popos=append(popos,patotes)
              papas=append(papas,rep(PiesFin[piss,2],length(patotes)))
            }
            
            papos=c()
            pospos=unique(popos)
            for(pipi in pospos){
              papos=append(papos,paste(papas[which(pipi == popos)], collapse=";"))  
            }
            
            write(paste(PasteApe(OriSeqNameIn,seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
            
            ###TODO: ADD info IN HTML here
            ###Test NewSequenceViewer = function(title,div_id,sequence,patterns,colors,tooltips,dflegends,starseq="",endseq=""){
            #newdf=data.frame(Type=c("Start","End","Multiple annotations"),Color=c("green","red","grey"))
            ##newdf=data.frame(Seqs=patseqs,Description=(regions)[,3],Color=(regions)[,4])
            ##newdf=unique(newdf)
            legdf=data.frame(Type=c("piRNA site", "RE site"),Color=c("orange","purple"))

            
            HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),NewSequenceViewer("","oldestseq",seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaI","SapI","SapI","Esp3I","Esp3I"),legdf,"",""),"<br>")
          }else{
            write(paste(PasteApe(OriSeqNameIn,seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
            legdf=data.frame(Type=c("RE site"),Color=c("purple"))
            ###TODO: ADD info IN HTML here
            HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),NewSequenceViewer("","oldestseq",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaI","SapI","SapI","Esp3I","Esp3I"),legdf,"",""),"<br>")
          }
        })
        
        ##ApeButtonDownOriginal
        output$downloadApeOriseq <- renderUI({
          downloadButton('DownOriApe', 'Download annotations on original sequence (Genbank)')
        })
        
        ##PiTab
        output$OriPiTab <- renderTable({
          piesinseq=Strfindpies(toupper(seqiqi),Pies,4)
          #Subsampling sequence and creating table
          pimattab=c()
          for(pipat in piesinseq){
            
            #matchpattern to find targets
            matseqpie=as.character(matchPattern(DNAString(pipat),DNAString(toupper(seqiqi)),max.mismatch=4,fixed=T))
            
            for(mat in matseqpie){
              pimattab=rbind(pimattab,cbind(pipat,mat,stringdist(pipat,mat,"hamming")))
            }
          }
          if(length(piesinseq) > 0){
            pimattab=cbind(PiesFin[pimattab[,1],2],pimattab)
            pimattab[,2]=paste(pimattab[,2],"|<-|",sep="")
            rownames(pimattab)=1:nrow(pimattab)
            colnames(pimattab)=c("piRNA locus","21-U reverse complement sequence","Matching sequence","Edit distance")
          }
          pimattab
        })
        
      }
      output$button4later <- renderUI({actionButton("actionRESET", label = "RESET")})
      ErrorFlag = 2
    }
    
    
    ###TWIST MODE
    ##Twist mode supercedes everything else
    if(ErrorFlag == 5){
      ##Load all results in dynamic ui
      output$DynamicUserInterface <- renderUI({uiOutput("AllResults")})
      
      
      output$AllResults <- renderUI({
        fluidRow(
          #actionButton("actionRESET", label = "RESET"),
          uiOutput("button4later"),
          htmlOutput("oldsequence"),
          #uiOutput("downloadApeOriseq"),
          hr()
          #dataTableOutput("OriPiTab")
        )
      })
      
      #output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
      output$button4later <- renderUI({HTML("<b></b>")})
      
      ##Perform analysis
      TwistResults=TwistSynthesisWithCoordinates(toupper(secdnausr))
      Twistcheck=TwistResults[[1]]
      
      if(Twistcheck == "pass"){
        output$oldsequence <-renderUI({
        HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br><b>Simple synthesis</b><br>",NewSequenceViewer("","oldestseq",secdnausr,c(),c(),c(),c(),"",""))
        })
        }else{
          
          errors=TwistResults[[2]]
          regions=TwistResults[[3]]
          
          starts=as.integer(regions[,1])
          ends=as.integer(regions[,2])
          
          patseqs=as.character(extractAt(DNAString(toupper(secdnausr)),IRanges(start=starts,end=ends)))
          newdf=data.frame(Seqs=patseqs,Description=(regions)[,3],Color=(regions)[,4])
          
          newdf=unique(newdf)
          
          legdf=data.frame(Type=as.character(newdf[,2]),Color=as.character(newdf[,3]))
          ###If original sequence
        output$oldsequence <-renderUI({
         #write(paste(PasteApe("Original_sequence",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
            
          HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br><b>Complex synthesis:</b><br><h5>",paste0(paste0(errors),sep="<br>",collapse=""),"<br></h5>",NewSequenceViewer("","oldestseq",secdnausr,c(as.character(newdf[,1])),c(as.character(newdf[,3])),c(as.character(newdf[,2])),legdf,"",""))
        })
        
        ##ApeButtonDownOriginal
        #output$downloadApeOriseq <- renderUI({
        #  downloadButton('DownOriApe', 'Download annotations on original sequence (Genbank)')
        #})
        
        
      }
      output$button4later <- renderUI({actionButton("actionRESET", label = "RESET")})
    }
    ####Where the future lays, i.e., the place where the sequence is created
    ####################################From here change the dynamic output so a table or something appears on the meanwhile, you can reander the gif for now
    
    if(ErrorFlag == 0){ ##Main routine
      ###Here's where the dynamic interface starts
      ##Render a data table, displaying current execution time
      ##Or expected completion time and status of process
      ##Status
      #output$DynamicUserInterface <- renderUI({verbatimTextOutput("Counter")})
      
      output$DynamicUserInterface <- renderUI({
        #HTML("")
        includeHTML("www/clock.html")
        })
      
      #ElapsedTime <- c(0)
      #TimerActive <- c(TRUE)
      
      ###Send custom messages to start reading status
      session$sendCustomMessage("submitted-job", TRUE)
      
      myFuture <- future({
        #withProgress(message = 'Generating transgene', style = "notification", detail = "(~4 min per kb for piRNA optimization and GLO algorithm)", value = 0, {
        ###Internal parameters
        RetrieveTop=5000
        
        writeLines(text="Performing Ribosomal optimization",con=statfiletowork)
        
        ######1st. step: Ribosomal Optimization
        if(FlaRi){ ###Ribosomal binding optimization
          testSeq=paste(c(seqDNA[1:39]),sep="",collapse="")
          write(paste("AAAA",testSeq,sep="",collapse=""),paste("DATA/users/",session_id,"/gene.fasta", sep=""))
          for(h in 1:100){
            testSeq2=repcds(testSeq,CAIS,AAtoCodF,8)
            write(paste("AAAA",testSeq2,sep="",collapse=""),paste("DATA/users/",session_id,"/gene.fasta", sep=""),append=T)
          }
          system(paste("sh bin/rnafold.sh",paste("DATA/users/",session_id,sep=""),RetrieveTop))
          inseqs=readLines(paste("DATA/users/",session_id,"/seqswithoutaaaas.txt", sep=""))
          fwd=grep("GGTCTC",x=inseqs)
          rev=grep("GAGACC",x=inseqs)
          all=unique(c(fwd,rev))
          if((length(all) != 0)&(length(all) != length(inseqs))){inseqs=inseqs[-c(all)]}
          id=order(sapply(inseqs,function(x){Strcondepies(x,Pies,4)}))[1]
          SeqStart=inseqs[id]
        }else{ ###Not optimization
          if(CodonAl != 6){
            SeqStart=repcds(paste(c(seqDNA[1:39]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl)
          }
          if(FlaModOnly){
            SeqStart=paste(c(seqDNA[1:39]),sep="",collapse="")
            }
        }
        
        
        #incProgress(3/10)
        
        ######2nd. step: Sequence Optimization
        writeLines(text="Adapting sequence by specified algorithm or codon usage",con=statfiletowork)
        
        ##CHeck box that overrides sequence optimization
        if(FlaModOnly){
          SeqtoOpt=paste(c(SeqStart,seqDNA[40:length(seqDNA)]),sep="",collapse="")
          }else{
        
        
        if((CodonAl <= 5)&(!(FlaPi))){ ###Use CAI equal to one
          SeqEnd=repcds(paste(c(seqDNA[40:length(seqDNA)]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl)
          
        }else{ 
          
          ##Check for Dans algoritmh, this should sepercede ribosomal optimization
          if(CodonAl == 6){ ###Produce start seq as originally it wasnt 
            
            write(paste(seqDNA,sep="",collapse=""),paste("DATA/users/",session_id,"/gene.fasta", sep=""))
            system(paste("perl bin/GLO_CLI_one_line.pl",paste("DATA/users/",session_id,"/gene.fasta", sep=""),">", paste("DATA/users/",session_id,"/GLO.fasta", sep="")))
            
            GLOseq=readLines(paste("DATA/users/",session_id,"/GLO.fasta", sep=""))
            GLOseq=unlist(strsplit(toupper(GLOseq),""))
            
            SeqStart=paste(c(GLOseq[1:39]),sep="",collapse="")
            SeqEnd=paste(c(GLOseq[40:length(GLOseq)]),sep="",collapse="")
            
          }else{ ###Use Codon Al as column of frequencies to mimic, and sample 100 times to reduce piRNAs   
            setrep=c()
            for(j in 1:100){
              setrep=append(setrep,repcds(paste(c(seqDNA[40:length(seqDNA)]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl))
            }
            inseqs=setrep
            fwd=grep("GGTCTC",x=inseqs)
            rev=grep("GAGACC",x=inseqs)
            all=unique(c(fwd,rev))
            if((length(all) != 0)&(length(all) != length(inseqs))){inseqs=inseqs[-c(all)]}
            id=order(sapply(inseqs,function(x){Strcondepies(x,Pies,4)}))[1]
            SeqEnd=inseqs[id]
          }
        }
        
        
        SeqtoOpt=paste(c(SeqStart,SeqEnd),sep="",collapse="")
        #############
        }
        #incProgress(3/10)
        
        ######3rd. step: Sequence maniputalion for piRNA removal
        writeLines(text="Removing piRNA sites",con=statfiletowork)
        
        ##If PiRNA removal
        if(FlaPi){
          pipipis=Strfindpies(SeqtoOpt,Pies,4)
          if(length(pipipis) > 0 ){
            stpos=c()
            for(pipi in pipipis){
              stpos=append(stpos,start(matchPattern(DNAString(pipi),DNAString(SeqtoOpt),max.mismatch=4,fixed=T)))
            }
            stpos=unique(c(stpos,stpos+3,stpos+6,stpos+9,stpos+12,stpos+15))
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,8)
            
            pipipis=Strfindpies(SeqtoOpt,Pies,4)
            if(length(pipipis) > 0 ){
              Iter=1
              nflag=TRUE
              while((Iter < 100)&(nflag)){
                
                stpos=c()
                for(pipi in pipipis){
                  stpos=append(stpos,start(matchPattern(DNAString(pipi),DNAString(SeqtoOpt),max.mismatch=4,fixed=T)))
                }
                stpos=unique(c(stpos,stpos+3,stpos+6,stpos+9,stpos+12,stpos+15))
                SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,8)
                
                pipipis=Strfindpies(SeqtoOpt,Pies,4)
                
                if(length(pipipis) > 0 ){Iter=1+Iter}else{nflag=FALSE}
              }
              ##Error Message for iterations
              if(Iter==100){output$ErrorMessage <- renderText({paste("Error: PiRNA removal did not work even after 100 iterations. Final number of piRNA sites found was: ",length(pipipis))})}
            }
          }
          
        }
        
        #incProgress(3/10)
        
        ######4th. step: Sequence maniputalion for restriction site removal
        writeLines(text="Restriction Sites removal",con=statfiletowork)
        
        Inenzy=c(as.character(gegegenzymes),as.character(gogogonzymes))
        
        #If restriction sites
        if(FlaEnz & (length(Inenzy) > 0)){
          enpat=c(as.character(enzy[Inenzy,"Site"]))
          for(pattemp in enpat){
            enpat=append(enpat,as.character(reverseComplement(DNAString(as.character(pattemp)))))
          }
          enpat=unique(enpat)
          #First search
          all=c()
          for(patito in enpat){
            all=append(all,grep(as.character(patito),x=SeqtoOpt))
          }
          all=unique(all)
          if(length(all) > 0 ){ ##Do proper biostrings match
            stpos=c()
            for(patito in enpat){
              stpos=append(stpos,start(matchPattern(DNAString(as.character(patito)),DNAString(SeqtoOpt),fixed=F)))
            }
            stpos=unique(c(stpos))
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,8)
            
            all=c()
            for(patito in enpat){
              all=append(all,grep(as.character(patito),x=SeqtoOpt))
            }
            all=unique(all)
            
            if(length(all) > 0 ){
              Iter=1
              nflag=TRUE
              while((Iter < 100)&(nflag)){
                stpos=c()
                for(patito in enpat){
                  stpos=append(stpos,start(matchPattern(DNAString(as.character(patito)),DNAString(SeqtoOpt),fixed=F)))
                }
                stpos=unique(c(stpos))
                SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,8)
                
                pipipis=Strfindpies(SeqtoOpt,Pies,4)
                
                all=c()
                for(patito in enpat){
                  all=append(all,grep(as.character(patito),x=SeqtoOpt))
                }
                all=unique(all)
                
                if(length(all) > 0 ){Iter=1+Iter}else{nflag=FALSE}
              }
              ##Error Message for iterations
              if(Iter==100){output$ErrorMessage <- renderText({paste("Error: Restriction site removal did not work even after 100 iterations. Final number of sites found was: ",length(all))})}
            }
          }
          
          
        }
        
        #incProgress(1/20)
        
        ######5th. step: Sequence maniputalion to add introns
        writeLines(text="Adding introns",con=statfiletowork)
        
        finalvec=c()
        #If introns
        if(FlaIn){
          
          finalvec=unlist(strsplit(toupper(SeqtoOpt),""))
          stpos=c()
          inpos=c()
          if(typedistint == 1){
            if(Flaframeint){
              stpos=append(stpos,start(matchPattern(DNAString("AGR"),DNAString(SeqtoOpt),fixed=F)))
              stpos=stpos + 1
              if(sum((stpos %% 3)==0) > 3){stpos=stpos[(stpos %% 3)==0]}else{
                stpos=append(stpos,start(matchPattern(DNAString("GR"),DNAString(SeqtoOpt),fixed=F)))
                stpos=stpos[(stpos %% 3)==0]
              }
              if(length(stpos)>3){
                if(sum((stpos > 50)&(stpos<150))>1){inpos=append(inpos,c(stpos[(stpos > 50)&(stpos<150)])[1])}else{inpos=c(51)}
                if(sum(stpos > (inpos[1]+150))>1){inpos=append(inpos,c(stpos[(stpos > (inpos[1]+150))])[1])}else{inpos=c(inpos[1],inpos[1]+150)}
                if(sum(stpos > (inpos[2]+150))>1){inpos=append(inpos,c(stpos[(stpos > (inpos[2]+150))])[1])}else{inpos=append(inpos,sample((inpos[2]+60):(length(finalvec)-1),1))}
              }else{
                inpos=c(99,249,399)
              }
              inposis=inpos[order(inpos)]
            }else{
              stpos=append(stpos,start(matchPattern(DNAString("AGR"),DNAString(SeqtoOpt),fixed=F)))
              stpos=stpos + 1
              if(length(stpos)>3){
                if(sum((stpos > 50)&(stpos<150))>1){inpos=append(inpos,c(stpos[(stpos > 50)&(stpos<150)])[1])}else{inpos=c(50)}
                if(sum(stpos > (inpos[1]+150))>1){inpos=append(inpos,c(stpos[(stpos > (inpos[1]+150))])[1])}else{inpos=c(inpos[1],inpos[1]+150)}
                if(sum(stpos > (inpos[2]+150))>1){inpos=append(inpos,c(stpos[(stpos > (inpos[2]+150))])[1])}else{inpos=append(inpos,sample((inpos[2]+50):(length(finalvec)-1),1))}
              }else{
                inpos=c(100,250,400)
              }
              inposis=inpos[order(inpos)]
            }
          }else{
            if(Flaframeint){
              stpos=append(stpos,start(matchPattern(DNAString("AGR"),DNAString(SeqtoOpt),fixed=F)))
              stpos=stpos + 1
              if(sum((stpos %% 3)==0) > 5){stpos=stpos[(stpos %% 3)==0]}else{
                stpos=append(stpos,start(matchPattern(DNAString("GR"),DNAString(SeqtoOpt),fixed=F)))
                stpos=stpos[(stpos %% 3)==0]
              }
              inpos=quantile(stpos,names=F)[c(2,3,4)]
              inposis=inpos[order(inpos)]
            }else{
              stpos=append(stpos,start(matchPattern(DNAString("AGR"),DNAString(SeqtoOpt),fixed=F)))
              stpos=stpos + 1
              inpos=quantile(stpos,names=F)[c(2,3,4)]
              inposis=inpos[order(inpos)]
            }
          }
          
          SeqtoOpt=paste(c(finalvec[1:inposis[1]],as.character(IntronSeqs[typeIn,1]),finalvec[(inposis[1]+1):inposis[2]],as.character(IntronSeqs[typeIn,2]),finalvec[(inposis[2]+1):inposis[3]], as.character(IntronSeqs[typeIn,3]),finalvec[(inposis[3]+1):length(finalvec)]),sep="",collapse="")
        }
        
        ########
        ##6th Addition of UTRs
        ######6th. step: Sequence maniputalion to add UTR
        #writeLines(text="Adding UTRs",con=statfiletowork)
        ##3Better toi add at the end
        #If UTRs
        ##if(FlaTURS){
        #  SeqtoOpt=paste(befUTR,SeqtoOpt,aftUTR,sep="",collapse="")
        #}
        
        
        ########Mentiras mias, segun yo identifique en que momento terminaba el futuro
        #incProgress(1/20)
        
        results=list(SeqtoOpt,finalvec)
        ####################################################################################################################
        results
        ###Consider writing SeqtoOPt in usr folder so no race is invoked
      }) 
      
      
      then(myFuture, onFulfilled = function(resultslist)
      
      {#############################################Future happening############################################################3
        #############################################################
        ######6th. step: Display of results
        SeqtoOpt=unlist(resultslist[1])
        finalvec=unlist(resultslist[2])
        
        ##Variables to consider
        #IF: FlaTURS
        #  befUTR,SeqtoOpt,aftUTR
        
        ##If annotation, show results
        #SeqtoOpt <- .
        ###Also, it seems future forgets about all variables created previously which makes sense
        ###This allows to add dynamic outputs
        output$DynamicUserInterface <- renderUI({uiOutput("AllResults")})
        ###Send custom messages to stop status
        session$sendCustomMessage("submitted-job", FALSE)
        
        ##Extra things coming from before
        gegegenzymes=input$Genzymes
        gogogonzymes=input$Oenzymes
        
        ###Clean status file
        if(file.exists(statfiletowork)){file.remove(statfiletowork)}
        
        Inenzy=c(as.character(gegegenzymes),as.character(gogogonzymes))
        
        ##This creates these dynamic outputs
        ##A clear function will maybe introduce errors asthe dynamic output differs, similarly, if those are clear just before getting displayed, I dont think will make a big difference
        
        
        ###FlaAnno will be always true, so clear function is fine, let's see now if clear function can work before optimization
        
        if(FlaAno){
          output$AllResults <- renderUI({
            fluidRow(
              #actionButton("actionRESET", label = "RESET"),
              #actionButton("actionRESET", label = "RESET"),
              uiOutput("button4later"),
              hr(),
              ##First optimized
              #DT::dataTableOutput('OriPiTab'),
              htmlOutput("newsequence"),
              uiOutput("downloadApeOptiseq"),
              uiOutput("downloadoptseq"),
              tableOutput("OptiPiTab"),
              ##Then previous
              htmlOutput("oldsequence"),
              uiOutput("downloadApeOriseq"),
              tableOutput("OriPiTab")
              #DT::dataTableOutput('OptiPiTab')
            )
          })
        }
        
        #output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
        output$button4later <- renderUI({HTML("<b></b>")})
        output$ErrorMessage <- renderText({""})
        
    
        ##Now sequence manipulation
        aaaads=""
        if((FlaRi)|(Fla5paaa)){aaaads="aaaa"}
        
        if(!(FlaIn)){finalvec=unlist(strsplit(toupper(SeqtoOpt),""))}
        if(!(FlaIn)){optsec=SeqtoOpt}else{optsec=paste(finalvec,sep="",collapse="")}
        
        output$downloadoptseq <- renderUI({
          if((ErrorFlag == 0) & !is.null(SeqtoOpt)) {
            if(CodonAl == 6){
              optsin="GLO"
            }else{
              optsin=colnames(CAIS)[CodonAl]
            }
            if(FlaRi){optsin=paste(optsin,"_OptimalRibosomalBinding",sep="",collapse="")}
            if(FlaPi){optsin=paste(optsin,"_RemovepiRNAHomology",sep="",collapse="")}
            if(FlaEnz & (length(Inenzy) > 0)){
              noenz=paste("_no",gsub("\\.","",Inenzy),sep="")
              optsin=paste(optsin,noenz,sep="",collapse="")
            }
            if(FlaIn){optsin=paste(optsin,"_withSyntheticIntrons",sep="",collapse="")}
            if(FlaTURS){optsin=paste(optsin,"_andUTRs",sep="",collapse="")}
            #IF: FlaTURS
            #  befUTR,SeqtoOpt,aftUTR
            
            #Write optimized sequence as a gene block
            write(paste(paste("-User_",unlist(strsplit(as.character(Sys.time()), " "))[2],sep=""),"CDS","gold",SeqtoOpt,sep="\t"),paste("DATA/users/",session_id,"/UserElements.tsv", sep=""), append=T)
            #Write Optimized sequence to fasta
            
            ###TODO: Change name in the fasta sequence, though realize that optsin will need to be changed when the codons usages are modified/removed
            write(paste(">Optimized_",SeqNameIn,":Codon-",optsin,"\n",befUTR,aaaads,SeqtoOpt,aftUTR,"\n",sep="",collapse=""),paste("DATA/users/",session_id,"/SeqOpop.fasta", sep=""))
            downloadButton('DownSeqOut', 'Download optimized sequence (Fasta)')
          }
        })
        
        SeqtoOpt=paste(befUTR,aaaads,SeqtoOpt,aftUTR,sep="",collapse="")
        ###Graphical Output
        ##Check for FlaPi as we would use that flag to add the piRNAs or not
        if(FlaAno){
          if((ErrorFlag == 0) & !is.null(SeqtoOpt)) {
            if(as.integer(typinput) == 1){ ###If original sequence
              output$oldsequence <-renderUI({
                seqiqi=toupper(seqiqi)
                
                piss=Strfindpies(seqiqi,Pies,4)
                if((length(piss)>0)&(FlaPi)){
                  popos=c()
                  papas=c()
                  for(pipi in piss){
                    patotes=as.character(matchPattern(DNAString(pipi),DNAString(seqiqi),max.mismatch=4,fixed=T))
                    popos=append(popos,patotes)
                    papas=append(papas,rep(PiesFin[piss,2],length(patotes)))
                  }
                  
                  papos=c()
                  pospos=unique(popos)
                  for(pipi in pospos){
                    papos=append(papos,paste(papas[which(pipi == popos)], collapse=";"))  
                  }
                  ##Write annotated file
                  write(paste(PasteApe(OriSeqNameIn,seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
                  
                  ###TODO: ADD info IN HTML here
                  #"<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),SequenceViewer(""
                  legdf=data.frame(Type=c("piRNA site", "RE site"),Color=c("orange","purple"))
                  
                  
                  HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),NewSequenceViewer("","oldestseq",seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaI","SapI","SapI","Esp3I","Esp3I"),legdf,"",""),"<br>")
                  
                  #HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),SequenceViewer("","oldestseq",seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaI","SapI","SapI","Esp3I","Esp3I")))
                }else{
                  write(paste(PasteApe(OriSeqNameIn,seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
                  
                  ###TODO: ADD info IN HTML here
                  #HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),SequenceViewer("","oldestseq",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaI","SapI","SapI","Esp3I","Esp3I")))
                  legdf=data.frame(Type=c("RE site"),Color=c("purple"))
                  ###TODO: ADD info IN HTML here
                  HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),NewSequenceViewer("","oldestseq",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaI","SapI","SapI","Esp3I","Esp3I"),legdf,"",""),"<br>")
                  
                  }
              })
              
              
              
              ##ApeButtonDownOriginal
              output$downloadApeOriseq <- renderUI({
                downloadButton('DownOriApe', 'Download annotations on original sequence (Genbank)')
              })
              
              ##PiTab
              if(FlaPi){
              output$OriPiTab <- renderTable({
                piesinseq=Strfindpies(toupper(seqiqi),Pies,4)
                #Subsampling sequence and creating table
                pimattab=c()
                for(pipat in piesinseq){
                  
                  #matchpattern to find targets
                  matseqpie=as.character(matchPattern(DNAString(pipat),DNAString(toupper(seqiqi)),max.mismatch=4,fixed=T))
                  
                  for(mat in matseqpie){
                    pimattab=rbind(pimattab,cbind(pipat,mat,stringdist(pipat,mat,"hamming")))
                  }
                }
                if(length(piesinseq) > 0){
                  pimattab=cbind(PiesFin[pimattab[,1],2],pimattab)
                  pimattab[,2]=paste(pimattab[,2],"|<-|",sep="")
                  rownames(pimattab)=1:nrow(pimattab)
                  colnames(pimattab)=c("piRNA locus","21-U reverse complement sequence","Matching sequence","Edit distance")
                }
                pimattab
              })
              }
              
            }
            
            ###For output sequence
            output$newsequence <-renderUI({
              coolpatterns=c()
              stpos=c()
              edpos=c()
              enpat=c()
              enpat=append(enpat,c("GGTCTC","GAGACC"))
              enpat=append(enpat,c("GCTCTTC","GAAGAGC"))
              enpat=append(enpat,c("CGTCTC","GAGACG"))
              colocolo=rainbow(length(enpat))
              
              for(t in 1:length(enpat)){
                stpos=start(matchPattern(DNAString(as.character(enpat[t])),DNAString(paste(finalvec,sep="",collapse="")),fixed=T))
                edpos=end(matchPattern(DNAString(as.character(enpat[t])),DNAString(paste(finalvec,sep="",collapse="")),fixed=T))
                if(length(stpos)>0){
                  coolpatterns=rbind(coolpatterns, cbind(stpos,edpos,rep(colocolo[t],length(stpos)),rep(enpat[t],length(stpos))))
                }
              }
              
              piss=Strfindpies(SeqtoOpt,Pies,4)
              if((length(piss)>0)&(FlaPi)){
                popos=c()
                papas=c()
                for(pipi in piss){
                  patotes=as.character(matchPattern(DNAString(pipi),DNAString(SeqtoOpt),max.mismatch=4,fixed=T))
                  popos=append(popos,patotes)
                  papas=append(papas,rep(PiesFin[piss,2],length(patotes)))
                }
                
                papos=c()
                pospos=unique(popos)
                for(pipi in pospos){
                  papos=append(papos,paste(papas[which(pipi == popos)], collapse=";"))  
                }
                
                ##Write annotated file
                write(paste(PasteApe(paste("Optimized_",OriSeqNameIn,sep=""),SeqtoOpt,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqpop.gb", sep=""))
                
                
                ###TODO: ADD info IN HTML here
                legdf=data.frame(Type=c("piRNA site", "RE site"),Color=c("orange","purple"))
                #"<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),SequenceViewer(""
                HTML("<br><b><h3>Optimized",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(optsec,CAIS,5,AAtoCodF),NewSequenceViewer("","newseq",SeqtoOpt,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaI","SapI","SapI","Esp3I","Esp3I"),legdf,"",""),"<br>")
                #legdf=data.frame(Type=c("piRNA site", "RE site"),Color=c("orange","purple"))
                
                
                #HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),NewSequenceViewer("","oldestseq",seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaI","SapI","SapI","Esp3I","Esp3I"),legdf,"",""),"<br>")
                
                
                }else{
                write(paste(PasteApe(paste("Optimized_",OriSeqNameIn,sep=""),SeqtoOpt,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqpop.gb", sep=""))
                  legdf=data.frame(Type=c("RE site"),Color=c("purple"))
                ###TODO: ADD info IN HTML here
                HTML("<br><b><h3>Optimized",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(optsec,CAIS,5,AAtoCodF),NewSequenceViewer("","newseq",SeqtoOpt,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaI","SapI","SapI","Esp3I","Esp3I"),legdf,"",""),"<br>")
                  #legdf=data.frame(Type=c("RE site"),Color=c("purple"))
                  ###TODO: ADD info IN HTML here
                  #HTML("<br><b><h3>",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),NewSequenceViewer("","oldestseq",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaI","SapI","SapI","Esp3I","Esp3I"),legdf,"",""),"<br>")
                  
                  }
            })
            
            output$downloadApeOptiseq <- renderUI({
              downloadButton('DownOptiApe', 'Download annotations on optimized sequence (Genbank)')
            })
            
            ##PiTab
            if(FlaPi){
            output$OptiPiTab <- renderTable({
              piesinseq=Strfindpies(toupper(optsec),Pies,4)
              #Subsampling sequence and creating table
              pimattab=c()
              for(pipat in piesinseq){
                
                #matchpattern to find targets
                matseqpie=as.character(matchPattern(DNAString(pipat),DNAString(toupper(optsec)),max.mismatch=4,fixed=T))
                
                for(mat in matseqpie){
                  pimattab=rbind(pimattab,cbind(pipat,mat,stringdist(pipat,mat,"hamming")))
                }
              }
              
              if(length(piesinseq) > 0){
                pimattab=cbind(PiesFin[pimattab[,1],2],pimattab)
                pimattab[,2]=paste(pimattab[,2],"|<-|",sep="")
                rownames(pimattab)=1:nrow(pimattab)
                colnames(pimattab)=c("piRNA locus","21-U reverse complement sequence","Matching sequence","Edit distance")
              }
              pimattab
            })
            }
          }
          
          output$button4later <- renderUI({actionButton("actionRESET", label = "RESET")})
          
        }### Perform annotation
      },
      onRejected = function(){
        inituiui()
        ###Send custom messages to stop status
        session$sendCustomMessage("submitted-job", FALSE)
        
        ###Clean status file
        if(file.exists(statfiletowork)){file.remove(statfiletowork)}
        
        output$ErrorMessage <- renderText({
          paste("Unxexpected error, please contact us with the details of your submitted job")
        })
        }
      ) ## end of main then
      
      ###End Main sequence adaptation routine
      ##Clean clock
      
      
      #}) ### end of progresses
      return(NULL)
    }
    
  })##ENd main observer 
  
  }) ####End of shiny server function
    
