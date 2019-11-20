#in Rstudio, ctrl+alt+enter send code to terminal pane
using DelimitedFiles
using RCall
using DataFrames
using Statistics
using Distributed
using ProgressMeter
addprocs(23)


tpmtableraw = readdlm("/home/users/sypark/00_Project/"*
            "01_thymoma/10_Final_data/"*
					  "01_expression/IRS4_corrected_v2/"*
					  "thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv");
tpmtableraw[1:5,1:5]
tpmtable=tpmtableraw[2:end,2:end];
tpmtablegene=tpmtableraw[:,1][2:end]
tpmtable=convert(AbstractArray{Float64,2},tpmtable)
logtpm = log10.(tpmtable.+1)
z_logtpm = similar(logtpm)
p_logtpm = similar(logtpm)

for i in 1:size(logtpm,2)
  z_logtpm[:,i] = (logtpm[:,i] .- mean(logtpm[:,i]))/std(logtpm[:,i])
  p_logtpm[:,i] = (logtpm[:,i])./sum(logtpm[:,i])
end

R"""
library(tidyverse);library(stringr)
gtf=read_tsv(paste0("/home/users/kjyi/Projects/thymus_single_cell/final/",
  "expression/IRS4/merged_gtf/Homo_sapiens.GRCh38.95.mod3.gtf"),skip=5,
  col_names = c(
    "chr","spircename","feature","start","end",
    "score","strand","frame","attribute"),
	col_types="ccciicccc") %>%
  filter(feature == "gene") %>%
	dplyr::select(chr,start,end,attribute) %>%
	mutate(attribute = str_replace(attribute, '.*gene_name \"','') %>%
	  str_replace('\".*','')) %>%
	as.data.frame()

chrsizes <- read_tsv('http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
					  col_names = c("chr","size"),col_types="ci") %>%
			 mutate(chr = str_replace(chr,'chrM','chrMT') %>% str_replace('chr','')) %>%
			 as.data.frame()
chrsizes <- chrsizes[1:24,] %>% arrange(as.numeric(chr))
""";
@rget gtf
@rget chrsizes

function createinterval2(chrsz; w=2*10^7,overlap=0.2)
	o = []
	for i in unique(chrsz.chr)
		chr=i
		st=1
		ed=1+w
		oo = [chr st ed]
		intv=w*overlap
		while chrsz.size[chrsz.chr.==i][1]>ed
			st +=intv
			ed +=intv
			oo=vcat(oo,[chr st ed])
		end
		push!(o,oo)
	end
	o
end


function segtogene(ch,s,E)
	gtf.attribute[(gtf.chr.==ch) .& (gtf.start .<E) .& (gtf.end .>s)]
end

function getcoordexp(myinterval,tpmtable)
	o=[]
	ong=[]
	oavg=[]
	#n=@distributed (+) for mi in myinterval
	@showprogress for mi in myinterval
		chr=mi[1,1]
		oo=zeros(size(mi,1),size(tpmtable,2))
		oong=zeros(size(mi,1),size(tpmtable,2))	
		for i in 1:size(mi,1)
			genes=segtogene(chr,mi[i,2],mi[i,3])
			for j in 1:size(tpmtable,1)
				@inbounds begin
					if tpmtablegene[j] in genes
						oo[i,:] = oo[i,:] + tpmtable[j,:]
						oong[i,:] = oong[i,:] .+ 1
					end
				end
			end
		end
		push!(o,oo)
		push!(ong,oong)
		push!(oavg,oo./oong)
		1
	end
	o,ong,oavg
end

myinterval=createinterval2(chrsizes,w=2*10^7,overlap=0.2);
myinterval[1]

function main(interval, table, prefix)
  o,ong,oavg=getcoordexp(interval,table);
  R"""
    o=$o
    ong=$ong
    oavg=$oavg
    myinterval=$myinterval
    casenames=$(tpmtableraw[1,2:end]) %>% unlist()
    chrnames=lapply(myinterval,function(x)x[1,1]) %>% unlist()
    names(o)=chrnames
    names(ong)=chrnames
    names(oavg)=chrnames
    for(i in 1:length(o)){
    	colnames(o[[i]]) = casenames
  	  colnames(ong[[i]]) = casenames
  	  colnames(oavg[[i]]) = casenames
  	  rownames(o[[i]]) = paste0(myinterval[[i]][,2],"_",myinterval[[i]][,3])
  	  rownames(ong[[i]]) = paste0(myinterval[[i]][,2],"_",myinterval[[i]][,3])
  	  rownames(oavg[[i]]) = paste0(myinterval[[i]][,2],"_",myinterval[[i]][,3])
    }
    write_rds(o,$(prefix*"_sum.Rds"))
    write_rds(ong,$(prefix*"_ngenes.Rds"))
    write_rds(oavg,$(prefix*"_average.Rds"))
  """;
  0
end

run(`mkdir -p data/expbycoord`)

main(myinterval, z_logtpm, "data/expbycoord/zlogtpm_w2x10^7_o_0.2")
main(myinterval, p_logtpm, "data/expbycoord/plogtpm_w2x10^7_o_0.2")
main(myinterval, tpmtable, "data/expbycoord/tpm_w2x10^7_o_0.2")
main(myinterval, logtpm, "data/expbycoord/logtpm_w2x10^7_o_0.2")
