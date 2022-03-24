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
GoldenE=enzy[which(enzy$GG =="TRUE"),]
OtherE=enzy[-c(which(enzy$GG=="TRUE")),]


##TO DO: 
#it will be better to make the data transformations and save them as a .RDS, for now let's keep them there

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
  \'search\': true,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
</script>                   
                   ")
    }else{
      patitos=matchPattern(DNAString(as.character("ATG")),DNAString(paste(sequence,sep="",collapse="")),fixed=T)[1]
      
      for(pat in patterns){
        patitos=c(patitos,matchPattern(DNAString(as.character(pat)),DNAString(paste(sequence,sep="",collapse="")),fixed=T))
      }
      
      
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
      
      subcol=c("green")
      subtol=c("ATG")
      
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
  \'search\': true,
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
  
  ############################################Make Ape file#################################
  PasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,tabibi,cai,list,PiesList){
    if (is.null(sequence)){return(NULL)}
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
    
    FileLines=append(FileLines,paste("COMMENT",paste("Codon Adaptation Index",as.character(CAIS)),sep="     "))
    FileLines=append(FileLines,paste("COMMENT",paste("%GC",as.character(GCp)),sep="     "))
    FileLines=append(FileLines,paste("COMMENT",paste("Total number of piRNAs targeting this sequence",as.character(NoPies)),sep="     "))
    
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
  
  ###########################################################################################
  ##############################Download handlers############################################
  ###########################################################################################
  
  ##################Download handler
  output$DownSeqOut <- downloadHandler(
    filename <- function() {
      paste("Genebuild", "fasta", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/SeqOpop.fasta", sep=""), file)
    },
  )
  
  output$DownBlockConstruct <- downloadHandler(
    filename <- function() {
      paste("Genebuild", "fasta", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/BlockConstruct.fasta", sep=""), file)
    },
  )
  
  
  output$DownOriApe <- downloadHandler(
    filename <- function() {
      paste("Input_sequence", "gb", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/Seqog.gb", sep=""), file)
    },
  )
  
  output$DownOptiApe <- downloadHandler(
    filename <- function() {
      paste("Output_sequence", "gb", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/Seqpop.gb", sep=""), file)
    },
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
        radioButtons("intypeinput", label = HTML("Input"),
                     choices = list("DNA" = 1, 
                                    "Protein" = 2), 
                     selected = 1, inline = TRUE, width='100%'),
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
                    choices = list("--- Codon table ---" = 1, "Ubiquitous" = 2, "Germline" = 3, "Neuronal" = 4, "Somatic" = 5, "--- Published algorithms ---" = 6, "Max. expression (Henrik Bringmann)" = 7, "Germline optimization (GLO, Dan Dickinson)" = 8), 
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
                        This option calculates and minimizes the folding energy of positions -4 to +39 which optimizes ribosomal binding. 
 Also, this option adds four 5' leading adenines (aaaaATG) to the optimized sequence. Note, this option will be turned off if Germline optimization is selected.
                    </div></p>"),
 fluidRow(
   column(8,
          checkboxInput("checkPirna", label = HTML("<b>Minimize <i>C. elegans</i> piRNA homology
                                                              [<a href=\"\" onclick=\"$('#explain_piRNA').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%')
   )),
 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_piRNA\">
                        This option minimizes sequence homology to piRNAs to improve germline expression (<a href=\"https://science.sciencemag.org/content/359/6375/587\">Zhang et al. (2018)</a>, <a href=\"https://academic.oup.com/nar/article/46/W1/W43/4979435\">Wu et al. (2018)</a>, <a href=\"http://cosbi4.ee.ncku.edu.tw/pirScan/\">pirScan</a>). Our algorithm removes sequences that have up to five mismatches with known piRNAs. This algorithm is computationally demanding and takes approximately four minutes per kilobase.
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
                                               c(OtherE$Enzyme), inline=TRUE),
                           checkboxInput("checkMetSites", label = HTML("<b>Dam/Dcm
                                                               [<a href=\"\" onclick=\"$('#explain_metsites').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
                           HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_metsites\">
                        Disable restriction sites affected by <a href=\"https://blog.addgene.org/plasmids-101-methylation-and-restriction-enzymes\">Dam/Dcm</a> methylases. </div></p>")
                    ),
                    #HTML("</tt>")
                  )),
 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_enzysites\">
                        Note, this optimization is performed after piRNA removal and may re-introduce piRNAs sites. Any piRNAs sites are annotated in the sequence output.
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
                                          "rps-20 (62 bp, 28% GC)" = 4
                                        ), 
                                        selected = 2, width='100%'),
                           radioButtons("intdistop",label = HTML("Intron placement"),
                                        choices = list("Early start" = 1, 
                                                       "Equi-distant" = 2), 
                                        selected = 1, width='100%', inline = TRUE),
                           checkboxInput("checkintframe", label = HTML("Force introns in reading frame"), value = FALSE)
          ))),
 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_introns\">
                        Introns can improve transgene expression (<a href=\"https://pubmed.ncbi.nlm.nih.gov/8244003/\">Okkema <i>et al.,</i> 1993</a>). The introns are indicated by lower-case letters and are inserted at consensus splice sites (<a href=\"https://www.ncbi.nlm.nih.gov/books/NBK20075/\"><i>C. elegans</i> II, 2. ed</a>). 
                        Introns within the first 150 basepairs are particularly efficient at improving germline expression (<a href=\"https://www.nature.com/articles/s41467-020-19898-0\">Al Johani <i>et al.,</i> 2020</a>).
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
                    </div></p>")),
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
    #verbatimTextOutput("PartialResult"),
    output$PartialResult <- renderText({})
    
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
    output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
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
  #####Trangene generation#####
  observeEvent(input$actionSeq, {  
    
    ErrorFlag=0
    
    CodonAl=as.integer(input$selectCAI)
    
    FlaPi=input$checkPirna
    FlaIn=input$checkIntron
    FlaRi=input$checkboxRibo
    FlaEnz=input$checkEnzySites
    FlaAna=input$checkAnal
    ###Patch in wrong way
    #FlaAno=input$checkAnno
    FlaAno=TRUE
    typinput=input$intypeinput
    secprousr=input$seqPROT
    secdnausr=input$seqDNA
    gegegenzymes=input$Genzymes
    gogogonzymes=input$Oenzymes
    
    ##Additional flag so inputs are obtain from the beginning of function
    if(FlaIn){
      typeIn=as.integer(input$intropt)
      typedistint=as.integer(input$intdistop)
      Flaframeint=input$checkintframe
    }
    
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
        ###Make sure to make analysis flag off
        FlaAna=FALSE
      }
    }else{
      seqDNA=unlist(strsplit(toupper(secdnausr),""))
    }
    
    ##Analytical mode supercedes codon table
    if(!(FlaAna)){
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
        paste("Warning: Sequence cannot remove piRNAs and mantain a Codon Adaptation Index equal to 1. piRNAs will not be removed")
      })
      FlaPi=FALSE
      updateCheckboxInput(session, "checkPirna", value = FALSE)
    } 
    
    if((ErrorFlag == 0) & ((FlaRi)&(CodonAl == 6))){ ##Check if GLO is required
      output$ErrorMessage <- renderText({
        paste("Warning: Sequence cannot be optimized for ribosome binding and Germline expression at the same time. GLO algoritmh will be run on full sequence")
      })
      FlaRi=FALSE
      updateCheckboxInput(session, "checkboxRibo", value = FALSE)
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
    
    ##Analytical mode supercedes main routine
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
      
      output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
      
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
            
            write(paste(PasteApe("Original_sequence",seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
            HTML(SequenceViewer("Original sequence","oldestseq",seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaI","SapI","SapI","Esp3I","Esp3I")))
          }else{
            write(paste(PasteApe("Original_sequence",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
            HTML(SequenceViewer("Original sequence","oldestseq",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaI","SapI","SapI","Esp3I","Esp3I")))
          }
        })
        
        ##ApeButtonDownOriginal
        output$downloadApeOriseq <- renderUI({
          downloadButton('DownOriApe', 'Download annotations on original sequence')
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
        }
        
        
        #incProgress(3/10)
        
        ######2nd. step: Sequence Optimization
        writeLines(text="Adapting sequence by specified algorithm or codon usage",con=statfiletowork)
        
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
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
            
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
                SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
                
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
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
            
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
                SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
                
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
        ##If annotation, show results
        #SeqtoOpt <- .
        ###Also, it seems future forgets about all variables created previously which makes sense
        ###THis allows to add dynamic outputs
        
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
        
        
        ###FlaAnno will be always true, so clear function is fine, let's see now if clear function can wrok before optimization
        
        if(FlaAno){
          output$AllResults <- renderUI({
            fluidRow(
              #actionButton("actionRESET", label = "RESET"),
              #actionButton("actionRESET", label = "RESET"),
              uiOutput("button4later"),
              verbatimTextOutput("PartialResult"),
              uiOutput("downloadoptseq"),
              hr(),
              htmlOutput("oldsequence"),
              uiOutput("downloadApeOriseq"),
              tableOutput("OriPiTab"),
              #DT::dataTableOutput('OriPiTab'),
              htmlOutput("newsequence"),
              uiOutput("downloadApeOptiseq"),
              tableOutput("OptiPiTab")
              #DT::dataTableOutput('OptiPiTab')
            )
          })
        }else{
          output$AllResults <- renderUI({
            fluidRow(
              #actionButton("actionRESET", label = "RESET"),
              uiOutput("button4later"),
              verbatimTextOutput("PartialResult"),
              uiOutput("downloadoptseq")
            )
          })
          
        }
        
        output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
        
        output$ErrorMessage <- renderText({""})
        aaaads=""
        if(FlaRi){aaaads="aaaa"}
        if(!(FlaIn)){optsec=SeqtoOpt}else{optsec=paste(finalvec,sep="",collapse="")}
        #PartialResult
        output$PartialResult <- renderText({
          paste(c("Results:\n","GC content ",CalculateGC(optsec),"%\n","CAI ",CalculateCAI(optsec,CAIS,5,AAtoCodF),"\nNot removed piRNA sites ",length(Strfindpies(optsec,Pies,4)),"\n",aaaads,SeqtoOpt),sep="",collapse="")
        })
        
        
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
            #Write optimized sequence as a gene block
            write(paste(paste("-User_",unlist(strsplit(as.character(Sys.time()), " "))[2],sep=""),"CDS","gold",SeqtoOpt,sep="\t"),paste("DATA/users/",session_id,"/UserElements.tsv", sep=""), append=T)
            #Write Optimized sequence to fasta
            write(paste(">Optimized_cDNA:Codon-",optsin,"\n",aaaads,SeqtoOpt,"\n",sep="",collapse=""),paste("DATA/users/",session_id,"/SeqOpop.fasta", sep=""))
            downloadButton('DownSeqOut', 'Download Optimized DNA sequence')
          }
        })
        
        ###Graphical Output
        if(FlaAno){
          if((ErrorFlag == 0) & !is.null(SeqtoOpt)) {
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
                  ##Write annotated file
                  write(paste(PasteApe("Original_sequence",seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
                  
                  HTML(SequenceViewer("Original sequence","oldestseq",seqiqi,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaI","SapI","SapI","Esp3I","Esp3I")))
                }else{
                  write(paste(PasteApe("Original_sequence",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
                  HTML(SequenceViewer("Original sequence","oldestseq",seqiqi,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaI","SapI","SapI","Esp3I","Esp3I")))
                }
              })
              
              
              
              ##ApeButtonDownOriginal
              output$downloadApeOriseq <- renderUI({
                downloadButton('DownOriApe', 'Download annotations on original sequence')
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
              
              
              if(!(FlaIn)){finalvec=unlist(strsplit(toupper(SeqtoOpt),""))}
              for(t in 1:length(enpat)){
                stpos=start(matchPattern(DNAString(as.character(enpat[t])),DNAString(paste(finalvec,sep="",collapse="")),fixed=T))
                edpos=end(matchPattern(DNAString(as.character(enpat[t])),DNAString(paste(finalvec,sep="",collapse="")),fixed=T))
                if(length(stpos)>0){
                  coolpatterns=rbind(coolpatterns, cbind(stpos,edpos,rep(colocolo[t],length(stpos)),rep(enpat[t],length(stpos))))
                }
              }
              
              piss=Strfindpies(SeqtoOpt,Pies,4)
              if(length(piss)>0){
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
                write(paste(PasteApe("Gene_contruct",SeqtoOpt,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies),collapse="\n"),paste("DATA/users/",session_id,"/Seqpop.gb", sep=""))
                
                
                
                HTML(SequenceViewer("Optimized sequence","newseq",SeqtoOpt,c(pospos,"GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c(rep("orange",length(pospos)),"purple","purple","purple","purple","purple","purple"),c(papos,"BsaI","BsaI","SapI","SapI","Esp3I","Esp3I")))
              }else{
                write(paste(PasteApe("Gene_contruct",SeqtoOpt,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaIrev","SapI","SapIrev","Esp3I","Esp3Irev"),CAIS,5,AAtoCodF,Pies),collapse="\n"),paste("DATA/users/",session_id,"/Seqpop.gb", sep=""))
                HTML(SequenceViewer("Optimized sequence","newseq",SeqtoOpt,c("GGTCTC","GAGACC","GCTCTTC","GAAGAGC","CGTCTC","GAGACG"),c("purple","purple","purple","purple","purple","purple"),c("BsaI","BsaI","SapI","SapI","Esp3I","Esp3I")))
              }
            })
            
            output$downloadApeOptiseq <- renderUI({
              downloadButton('DownOptiApe', 'Download annotations on optimized sequence')
            })
            
            ##PiTab
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
    
