library(data.table)
library(stringr)
dat = fread("./data/annotations/sample_annotations_net_validation.csv")
#some decorations were removed during data entry.  we fix here.
RNAi_name_fix = matrix(c("C15C74", "C15C7.4"   ,
"K02F25",   "K02F2.5" ,
"R1022",   "R102.2" ,
"Y44A6D2",  "Y44A6D.2", 
"Y9C2VA1", "Y9C2UA.1" ,
"casy1", "casy-1"   ,
"egl21","egl-21"    ,
"flp12",   "flp-12" ,
"flp14",    "flp-14" , 
"flp1","flp-1",
"flp5","flp-5",
"flp9","flp-9",
"ida1","ida-1",
"ins24","ins-24",
"ins26","ins-26",
"ins30","ins-30",
"ins6","ins-6",
"mec12","mec-12",
"nlp15","nlp-15",
"nlp21","nlp-21",
"nlp50","nlp-50",
"pdf1","pdf-1",
"pgal1","pgal-1",
"pghm1","pghm-1",
"sbt1","sbt-1",
"snet1","snet-1",
"snt4","snt-4",
"ttr29","ttr-29",
"zig2","zig-2"),ncol=2,byrow=T)

for ( i in 1:nrow(RNAi_name_fix)){
	dat[RNAi==RNAi_name_fix[i],]$RNAi = RNAi_name_fix[i,2]
	dat[IntendedRNAi==RNAi_name_fix[i],]$IntendedRNAi = RNAi_name_fix[i,2]
}





fwrite(dat, "./data/annotations/sample_annotations_net_validation_fixed.csv", quote = TRUE)
