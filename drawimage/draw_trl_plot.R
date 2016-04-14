library(rjson)
"%+%" <- function(a, b) paste0(a, b)
options(stringsAsFactors = FALSE)
setwd("~/workspace/r")
OUTDIR="/ifs/rcgroups/profile-scratch/xg015/projectOutputScartch/breakmerImg2/006xg2/res/plate"
setwd(OUTDIR)
JSON_F="/ifs/rcgroups/profile/xg015/workspace/breakmerImg2/006xg2/res/plate/Pos_15-N14208-FFPE_B3_L_015838_HC_000668_363/res/output/CCNE1/contig1.json"
CIGAR_STR="93M82S,82M"
OUT_F="/home/xg015/workspace/t.png"
#these functions are adapted from https://www.biostars.org/p/9335/
matcher <- function(pattern, x) {
  ind = gregexpr(pattern, x)[[1]]
  start = as.numeric(ind)
  end = start + attr(ind, "match.length")- 2
  apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2]));
}
doone <- function(c, cigar) {
  pat <- paste("\\d+", c , sep="")
  sum(as.numeric(matcher(pat, cigar)), na.rm=T)
}
cigarsums <- function(cigar, chars=c("M","N","D","I","S","H", "P", "X", "=")) {
  sapply (chars, doone, cigar)
}
#end of the code from biostars
setup_cigar_obj<-function(cigar_str){
  #cigar_str=CIGAR_STR
  cigar_str=strsplit(cigar_str, ",")[1]
  cigar_sums=cigarsums(cigar_str[[1]][1])
  cigar_sums=as.data.frame(t(as.data.frame(cigar_sums)))
  col_names=colnames(cigar_sums)
  #set up onlyOneMandS
  only_one_M_and_S=""
  t=as.data.frame(table(col_names))
  if(t[t$col_names=="M",]$Freq==1 & t[t$col_names=="S",]$Freq==1){only_one_M_and_S="T"}
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
  return (list(only_one_M_and_S=only_one_M_and_S, cigar_label=cigar_label))
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
draw<-function(grid_, cigar_obj, out_f){
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
  seekViewport("header")
  grid.text("genename", x=0.5, y=0.5, just = c("center", "bottom"), gp=gpar(fontsize=fontsize, col=color, family="century"), check=TRUE)
  #end of annotation
  dev.off()
  #end os draw -------------------------------------------------------------------------------------
}
cigar_obj=setup_cigar_obj(CIGAR_STR)
cigar_obj
grid_= setup_grid(JSON_F)
grid_
draw(grid_, cigar_obj, OUT_F)
