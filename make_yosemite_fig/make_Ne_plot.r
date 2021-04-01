# load libraries
library(snpReady)
library(vcfR)
library(ggplot2)

# list all files
files = list.files()
files = files[grep('vcf', files)]

# function that takes a string of the form 'm|n'
# and returns the numeric m+n
cast_as_numeric = function(str){
  return(sum(as.numeric(strsplit(str, '|')[[1]][c(1,3)])))
}

# function that takes a matrix of values of the form 'm|n'
# and returns a new matrix of their sums
convert_to_counts = function(str_mat){
  return(matrix(unlist(lapply(str_mat, cast_as_numeric)), ncol = ncol(str_mat)))
}

# function that takes a VCF filename and calculates the Ne for that file
calc_Ne = function(filename){
  # read file
  vcf = read.vcfR(filename)
  # get the genotype matrix
  mat = as.matrix(vcf@gt)
  # transform (so that individs are on rows, markers on cols
  mat = t(mat)
  # drop the first row (which indicates the formatting, e.g. 'GT')
  mat = mat[2:nrow(mat), ]
  # convert to a numeric matrix of allele counts
  cts = convert_to_counts(mat)
  colnames(cts) = seq(ncol(cts))
  # calculate and extract the Ne
  res = popgen(cts, subgroups=rep(1, nrow(cts)))
  Ne = res$whole$Variability[1,1]
  return(Ne)
} 

print(files)

ts = c()
Nes = c()
for (filename in files){
  t = as.numeric(strsplit(strsplit(filename, '_t-')[[1]][2], "_spp")[[1]][1])
  Ne = calc_Ne(filename)
  ts = c(ts, t)
  Nes = c(Nes, Ne)
}

print('ts')
print(ts)
print('Nes')
print(Nes)

df = data.frame(t=ts, Ne=Nes)
print(df)

Ne_plot = ggplot(data=df) +
  geom_line(aes(x=t, y=Ne))

ggsave(Ne_plot, file='Ne_plot.pdf',
       width=30, height=20, units='cm', dpi=1200)
