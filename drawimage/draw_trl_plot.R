args=commandArgs(trailingOnly=TRUE)
print (length(args))
if(length(args)!=9){
  stop("Usage: Rscript --vanilla draw_trl_plot.R reads_json_file cigar_string output_png gene1 gene1_position gene2 gene2_position contig_name aggregate_pair_name \nexample: Rscript --vanilla draw_trl_plot.R ./contig1.json 93M82S,82M ./example.png CCNE1 19:30308080 null 19:31061218 CCNE1_contig1 17442, Pos_15-N14208-FFPE_B3_L_015838_HC_000668_363-Normal_Liver_B3-N_L_015840_HC_000668_363")
}
if(0==length(args)){
#### test setting##############
#72242
JSON_F="/ifs/rcgroups/profile/xg015/workspace/breakmerImg2/006xg2/res/plate/Pos_15-N14208-FFPE_B3_L_015838_HC_000668_363/res/output/FANCD2/contig5.json"
CIGAR_STR="14S22M1I2D96M4D21M9D76M"
OUT_F="/ifs/rcgroups/profile/xg015/workspace/breakmerImg2/006xg2/res/plate/img/72242.png"
GENE1="FANCD2"
POSITION1="chr3:10089452"
GENE2="FANCD2"
POSITION2="chr3:10089452"
CONTIG="FANCD2_contig5"
AGNAME="Pos_15-N14208-FFPE_B3_L_015838_HC_000668_363-Normal_Liver_B3-N_L_015840_HC_000668_363"
####end of test setting##############
}else{
JSON_F=args[1]
CIGAR_STR=args[2]
OUT_F=args[3]
GENE1=args[4]
POSITION1=args[5]
GENE2=args[6]
POSITION2=args[7]
CONTIG=args[8]
AGNAME=args[9]
print(JSON_F)
print(CIGAR_STR)
print(OUT_F)
print(GENE1)
print(POSITION1)
print(GENE2)
print(POSITION2)
print(CONTIG)
print(AGNAME)
}
library(rjson)
library(grid)
"%+%" <- function(a, b) paste0(a, b)
options(stringsAsFactors = FALSE)
#deprecated
#these functions are adapted from https://www.biostars.org/p/9335/
matcher <- function(pattern, x) {
  ind = gregexpr(pattern, x)[[1]]
  start = as.numeric(ind)
  end = start + attr(ind, "match.length")- 2
  apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2]));
}
#deprecated
doone <- function(c, cigar) {
  pat <- paste("\\d+", c , sep="")
  sum(as.numeric(matcher(pat, cigar)), na.rm=T)
}
#deprecated
cigarsums <- function(cigar, chars=c("M","N","D","I","S","H", "P", "X", "=")) {
  sapply (chars, doone, cigar)
}
#end of the code from biostars
#deprecated
setup_cigar_obj_BK<-function(cigar_str){
  #cigar_str=CIGAR_STR
  cigar_str=strsplit(cigar_str, ",")[1]
  cigar_sums=cigarsums(cigar_str[[1]][1])
  cigar_sums=as.data.frame(t(as.data.frame(cigar_sums)))
  col_names=colnames(cigar_sums)
  mFirst=F
  foundS=F
  foundM=F
  for(t in 1:length(col_names)){
    print(t)
    if("M"==col_names[t]){foundM=T}
    if("S"==col_names[t]){foundS=T}
    if(foundM==T & foundS==F){
      mFirst=T
      break
    } 
    if(foundM==F & foundS==T){
      mFirst=F
      break
    } 
  }
  #set up onlyOneMandS
  only_one_M_and_S=F
  t=as.data.frame(table(col_names))
  if(t[t$col_names=="M",]$Freq==1 & t[t$col_names=="S",]$Freq==1){only_one_M_and_S=T}
  #setup cigar_label
  cigar_label=character()
  for(i in 1:length(col_names)){
    type=col_names[i]
    n=as.numeric(cigar_sums[i])
    if(n==0){next}
    for(j in 1:n){
      cigar_label=append(cigar_label,type, length(cigar_label))
    }    
  }
  res=list(only_one_M_and_S=only_one_M_and_S, cigar_label=cigar_label, mFirst=mFirst)
  return (res)
}
setup_cigar_obj<-function(cigar_str){
  cigar_str=CIGAR_STR
  cigar_count=unlist(strsplit(strsplit(cigar_str, ",")[[1]][1], "[^[:digit:]]"))
  cigar_name=unlist(strsplit(strsplit(cigar_str, ",")[[1]][1], "[[:digit:]]"))
  cigar_name=cigar_name[cigar_name!=""]
  mFirst=F
  foundS=F
  foundM=F
  for(t in 1:length(cigar_name)){
    if("M"==cigar_name[t]){foundM=T}
    if("S"==cigar_name[t]){foundS=T}
    if(foundM==T & foundS==F){
      mFirst=T
      break
    } 
  }
  #set up onlyOneMandS
  only_one_M_and_S=F
  if(table(cigar_name)["M"]==1){only_one_M_and_S=T}
  
  #setup cigar_label
  cigar_label=character()
  for(i in 1:length(cigar_name)){
    type=cigar_name[i]
    if(type=="D"){next}
    n=cigar_count[i]
    for(j in 1:n){
      cigar_label=append(cigar_label,type, length(cigar_label))
    }    
  }
  res=list(only_one_M_and_S=only_one_M_and_S, cigar_label=cigar_label, mFirst=mFirst)
  return (res)
}
setup_annotation_obj<-function(cigar_obj, gene1,position1,gene2,position2, contig){
  if(!cigar_obj$only_one_M_and_S){return (NULL)}
  gene1=GENE1
  position1=POSITION1
  gene2=GENE2
  position2=POSITION2
  contig=CONTIG
  
  mFirst=cigar_obj$mFirst  
  primarygene=""
  primarygeneposition=""
  partnergene=""
  partnergeneposition=""
  if(length(unlist(strsplit(contig, "_")))>1){primarygene=unlist(strsplit(contig, "_"))[1]}
  if(primarygene==""){return (NULL)}  
  if(primarygene==gene1){
    primarygeneposition=position1  
    partnergene=gene2  
    partnergeneposition=position2
  }else if(primarygene==gene2){
    primarygeneposition=position2  
    partnergene=gene1  
    partnergeneposition=position1
  }else{return (NULL)}
  res=list(primarygene=primarygene,
       primarygeneposition=primarygeneposition,
       partnergene=partnergene,
       partnergeneposition=partnergeneposition,
       mFirst=mFirst)
  return (res)
}
setup_grid<-function(json_f){
  #json_f=JSON_F
  contig_json <- fromJSON(file=json_f)
  reads=contig_json$reads
  contig_id=contig_json$contigId
  nOfReads=length(reads)
  start_pos=-1
  for(i in 1:nOfReads){
    d=reads[[i]]
    seq=d$seq
    if(i==1){
      start_pos=d$start
    }
  }
  contig_seq=strsplit(contig_json$contigSeq, "")[[1]]
  if(0<start_pos){contig_seq=append(contig_seq, 1:start_pos, after = 0)}
  nOfCols=length(contig_seq)+start_pos
  nOfRows=nOfReads+1
  grid_ <- matrix(, nrow = nOfRows, ncol = nOfCols )
  grid_[1,]<-contig_seq
  for(row in 2:(nOfRows)){
    read_row=row-1
    read=reads[[read_row]]
    read_char=strsplit(read$seq, "")[[1]]
    if(0<read$start){read_char=append(read_char, 1:read$start, after = 0)}
    end_len=nOfCols-length(read_char)
    if(0<end_len){read_char=append(read_char, 1:end_len, after = length(read_char))}
    grid_[row,] <-read_char 
  }
  grid_[grid_!="A" & grid_!="G" & grid_!="T" & grid_!="C" & grid_!="N"]<-NA
  return (grid_)
}
draw<-function(grid_, cigar_obj, out_f, annotation_obj){
  only_one_M_and_S=cigar_obj$only_one_M_and_S
  cigar_label=cigar_obj$cigar_label
  #out_f=OUT_F
  filename=out_f
  fontsize=12
  grid_nofrows=dim(grid_)[1]
  grid_nofcols=dim(grid_)[2]
  grid_w=grid_nofcols*7
  grid_h=grid_nofrows*12
  #draw -----------------------------------------------------------
  png(filename=filename, width=grid_w, height=grid_h)
  grid.newpage()
  top.vp <- viewport(layout=grid.layout(2, 1, widths=unit(1, "null"), heights=unit(2, c("lines","null"))))
  header <- viewport(layout.pos.col = 1, layout.pos.row = 1, name = "header")
  box <- viewport(layout.pos.col = 1, layout.pos.row = 2, name = "box")
  splot <- vpTree(top.vp, vpList(header, box))
  pushViewport(splot)
  seekViewport("box")
  #starting point x y
  x=0
  y=1  
  yInx=y
  for (i in 1:grid_nofrows) {
    yInx=yInx-1/grid_nofrows
    xInx=x
    for(j in 1:grid_nofcols){      
      base=grid_[i,j]
      contig_ele=grid_[1,j]
      if(!is.na(base)){
        color="black"
        type=cigar_label[j]
        if(i==1){
          color="black"
        }else if(type=="M" & base==contig_ele){
          color="aquamarine4"
        }else if(type=="S" & base==contig_ele){
          color="orange"
        }else if(type=="I" & base==contig_ele){
          color="red"
        }else if(base!=contig_ele){
          color="black"
        }else{
          color="grey"
        }
        grid.text(base, x=xInx, y=yInx, just = c("left", "bottom"), 
                  gp=gpar(fontsize=fontsize, col=color, family="century"), check=TRUE)
      }
      xInx=xInx+1/grid_nofcols
    }
  }
  #annotation
  if(!is.null(annotation_obj)){
    seekViewport("header")
    if(annotation_obj$mFirst){
      grid.text(annotation_obj$primarygene, x=0.1, y=0.2, just = c("left", "bottom"), gp=gpar(fontsize=fontsize, col="aquamarine4", family="century"), check=TRUE)
      grid.text(annotation_obj$partnergene, x=0.7, y=0.2, just = c("left", "bottom"), gp=gpar(fontsize=fontsize, col="orange", family="century"), check=TRUE)
    }else{
      grid.text(annotation_obj$primarygene, x=0.8, y=0.2, just = c("left", "bottom"), gp=gpar(fontsize=fontsize, col="aquamarine4", family="century"), check=TRUE)
      grid.text(annotation_obj$partnergene, x=0.1, y=0.2, just = c("left", "bottom"), gp=gpar(fontsize=fontsize, col="orange", family="century"), check=TRUE)
    }
  }
  #end of annotation
  dev.off()
  #end of draw -------------------------------------------------------------------------------------
}
######################main entry##################
cigar_obj=setup_cigar_obj(CIGAR_STR)
grid_= setup_grid(JSON_F)
annotation_obj=setup_annotation_obj(cigar_obj, GENE1,POSITION1,GENE2,POSITION2, CONTIG)
draw(grid_, cigar_obj, OUT_F, annotation_obj)


#test
#example
JSON_F="/ifs/rcgroups/profile/xg015/workspace/breakmerImg2/006xg2/res/plate/Pos_15-N14208-FFPE_B3_L_015838_HC_000668_363/res/output/CCNE1/contig1.json"
CIGAR_STR="93M82S,82M"
GENE1="CCNE1"
POSITION1="19:30308080"
GENE2=""
POSITION2="19:31061218"
CONTIG="CCNE1_contig1"
OUT_F="/home/xg015/workspace/t.png"
AGNAME="Pos_15-N14208-FFPE_B3_L_015838_HC_000668_363"
