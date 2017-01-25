rm(list = ls())

## SUMMARY
# read snps
# calculate dist mx
# fit histogram

# Imports
library(cluster)
library(Biostrings)
library(dendextend)

# Parameterize script
fn = "./h3da_snp.fasta"
ancestral = "AGAATACATTTTTCCACTACCA"
cuts = 6

# Load file
raw = readLines(fn)
df = data.frame(name="FOO", snps="BAR", stringsAsFactors=FALSE)

## Define subroutines
numChr = function(char, s){
    # retunrs number of char in s
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
}

isDNA = function(s){
    # returns true if string is all A, C, G, T, N
    l = nchar(s)
    a = numChr("A", s)
    c = numChr("C", s)
    g = numChr("G", s)
    t = numChr("T", s)
    n = numChr("N", s)
    if (a + g + c + t + n == l){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Extract the names and snp strings
k = 1
for (i in 1:length(raw)){
    if (substr(raw[i],1,1)==">") {
        n = raw[i]
        s = raw[i+1]
        if (isDNA(s)){
            if (!(s %in% df$snps)){
                df = rbind(df, c(n, s))
                k = k + 1
            }
        }
    }
}
df = df[-1,] # remove the foobar row from the initialization step above

# Build the distance matrix of dist from org to the ancestral
dists = pairwiseAlignment(df$snps, ancestral, scoreOnly=TRUE)
dists = dists - min(dists)
dmx = dist(dists)

# Build the tree
hc = hclust(dmx)
hc_cut = cutree(hc, k = cuts) # cut the tree at <cuts> leafs
hcd = as.dendrogram(hc)
hcd_color = color_branches(hcd, k = cuts)

# Plot the tree
pdf("kalosa-kenyon_h3da.pdf")
plot(hcd_color, main="Flowering date gene h3da in 3000 rice cultivars",
        ylab = paste("Distance from ancestral sequence", ancestral),
        xlab = "Cultivar (h3da variant)",
        leaflab = "none"
    )

# Plot emperical distribution
plot(
        ecdf(dists),
        main="Emperical distribution for h3da distance from ancestor",
        xlab=paste("Distance from", ancestral),
        ylab="Emperical CDF"
    )

# Create histogram that is binrange=range_of_tree_cut and binheight=rel_frq
# Is there a gap?
# NOTE: range must be min > 0
bks = c() # this is the vec/list for the bin sizes in the histogram
eps = 0.00001
for (i in unique(hc_cut)){
    r = range(dists[hc_cut == i])
    bks = c(bks, r[1]-eps, r[2]+eps)
}
bks = sort(bks)
hist(dists, breaks = bks, main = "Histogram of genetic distances",
        xlab = paste("Distance from", ancestral),
        ylab = "Frequency"
    )
dev.off()
