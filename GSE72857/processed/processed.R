library(igraph)

in_silico_perturbation <- function(pset, rset, up, dn, strategy="v1", quantile=0.5){
  
  #pset, expression matrix to be preturbed, tfs in rows and samples in columns
  #rset, refererence expression matrix, tfs in rows and samples in columns
  #tf_up/dn, character vector, tfs to be up/down-regulated
  #strategy, perturbation strategies, specifically,
  #          v1, TF-wise perturbation, sample values for up/dn TFs from corresponding TF
  #          v2, TF-wise perturbation, sample values for up/dn TFs from all TF
  #          v3, TF-vector-wise perturbation, sample TF-vector
  #quantile, top x quantile for selecting tf values
  
  genes <- intersect(rownames(pset), rownames(rset))
  pset <- pset[genes, ]
  rset <- rset[genes, ]
  up <- intersect(genes, up)
  dn <- intersect(genes, dn)
  
  if (strategy == "v1"){
    for (x in up) pset[x, ] <- sample(rset[x, ][rset[x, ] >= quantile(rset[x, ], (1-quantile))], ncol(pset), replace=T)
    for (x in dn) pset[x, ] <- sample(rset[x, ][rset[x, ] <= quantile(rset[x, ], quantile)], ncol(pset), replace=T)}
  else if (strategy == "v2"){
    for (x in up) pset[x, ] <- sample(rset[rset >= quantile(rset, (1-quantile))], ncol(pset), replace=T)
    for (x in dn) pset[x, ] <- sample(rset[rset <= quantile(rset, quantile)], ncol(pset), replace=T)}
  else if (strategy == "v3"){
    pset[c(up, dn), ] <- rset[c(up, dn), sample(1:ncol(rset), ncol(pset), replace=T)]}
  return(pset)}

####################################################
##########    generate perturbed z data   ##########
####################################################

get_all_targets <- function(graph, from, direction, to=unique(V(graph)$name), depth=Inf){
  
  if (length(from) != length(direction)) stop("matched from and direction required")
  if (length(setdiff(direction, c(-1, 1))) > 0) stop("direction should be 1 (positive) or -1 (negative)")
  
  sgn <- structure(direction, names=from)
  to <- setdiff(intersect(unique(V(graph)$name), to), from)
  gdist <- distances(graph, from, to, mode="out", weights=NA)
  gdist <- gdist[, apply(gdist, 2, min) < depth]
  
  paths <- list()
  for (t in colnames(gdist)){
    f <- rownames(gdist)[gdist[, t] == min(gdist[, t])]
    path <- list()
    for (ff in f){
      p <- all_shortest_paths(graph, ff, t, mode="out", weights=NA)$res
      path[[ff]] <- lapply(p, function(pp) as_ids(pp))}
    path <- structure(do.call(c, path), names=NULL)
    path <- unlist(lapply(path, function(p, graph, sgn){
      idx <- c()
      for (i in 1:(length(p)-1)) idx <- c(idx, get.edge.ids(graph, c(p[i], p[i+1])))
      structure(sgn[p[1]], names=NULL) * Reduce("*", E(graph)$direction[idx])
    }, graph=graph, sgn=sgn))
    paths[[t]] <- sum(path)/length(path)}
  paths <- c(sgn, unlist(paths))
  return(paths)}

##################################################################
##########    get target tfs from the trrust database   ##########
##################################################################

annotation <- read.csv("./annotation.txt", header=F, stringsAsFactors=F)$V1
tfs <- read.table("./tfs.txt", quote="\"", comment.char="", stringsAsFactors=F)$V1
data_z <- read.table("./data_z.csv.gz", quote="\"", comment.char="", stringsAsFactors=F)
data_z <- structure(t(data_z), dimnames=list(tfs, annotation))

trrust <- read.delim("./trrust_rawdata.mouse.tsv", header=F, stringsAsFactors=F)
trrust <- trrust[trrust$V2 %in% trrust$V1, ]
trrust <- trrust[trrust$V3 == "Activation" | trrust$V3 == "Repression", ]
trrust <- graph_from_data_frame(data.frame(from=trrust$V1, to=trrust$V2,
                                           direction=c(-1, 1)[factor(trrust$V3, levels=c("Repression", "Activation"))],
                                           stringsAsFactors=F))
trrust <- simplify(trrust, remove.multiple=T, remove.loops=F, edge.attr.comb="first")

set <- get_all_targets(trrust, from=c("Gata1", "Spi1"), direction=c(-1, 1), depth=6)
dset1 <- in_silico_perturbation(pset=data_z[, colnames(data_z) == "CMP"], rset=data_z,
                                up=names(set[set > 0]), dn=names(set[set < 0]), strategy="v1", quantile=0.01)
colnames(dset1) <- rep("GMP_P", ncol(dset1))

set <- get_all_targets(trrust, from=c("Gata1", "Spi1"), direction=c(1, -1), depth=6)
dset2 <- in_silico_perturbation(pset=data_z[, colnames(data_z) == "CMP"], rset=data_z,
                                up=names(set[set > 0]), dn=names(set[set < 0]), strategy="v1", quantile=0.01)
colnames(dset2) <- rep("MEP_P", ncol(dset2))

perturb_z <- cbind(data_z, dset1, dset2)
gz <- gzfile("../perturb/z_perturb.csv.gz", "w")
write.table(t(perturb_z), gz, quote=F, row.names=F, col.names=F)
close(gz)
write.table(colnames(perturb_z), file="../perturb/annotation_perturb.txt", quote=F, row.names=F, col.names=F)

#####################################################
##########    prepare perturbation data    ##########
#####################################################

annotation <- read.csv("./annotation.txt", header=F, stringsAsFactors=F)$V1
tfs <- read.table("./tfs.txt", quote="\"", comment.char="", stringsAsFactors=F)$V1
data_z <- read.table("./data_z.csv.gz", quote="\"", comment.char="", stringsAsFactors=F)
data_z <- structure(t(data_z), dimnames=list(tfs, annotation))

trrust <- read.delim("./trrust_rawdata.mouse.tsv", header=F, stringsAsFactors=F)
trrust <- trrust[trrust$V2 %in% trrust$V1, ]
trrust <- trrust[trrust$V3 == "Activation" | trrust$V3 == "Repression", ]
trrust <- graph_from_data_frame(data.frame(from=trrust$V1, to=trrust$V2,
                                           direction=c(-1, 1)[factor(trrust$V3, levels=c("Repression", "Activation"))],
                                           stringsAsFactors=F))
trrust <- simplify(trrust, remove.multiple=T, remove.loops=F, edge.attr.comb="first")

p <- unlist(lapply(rownames(data_z), function(x, dat){
  wilcox.test(dat[x, colnames(dat) == "GMP"], dat[x, colnames(dat) == "MEP"])$p.value
}, dat=data_z))
d <- rowMeans(data_z[, colnames(data_z) == "GMP"]) - rowMeans(data_z[, colnames(data_z) == "MEP"])
names(p) <- names(d) <- rownames(data_z)
dn <- names(sort(d, decreasing=F))[1:5]
up <- names(sort(d, decreasing=T))[1:5]
sets <- apply(expand.grid(dn, up), 1, function(x, g) get_all_targets(g, from=x, direction=c(-1, 1), depth=6), g=trrust)
names(sets) <- apply(expand.grid(dn, up), 1, function(x) paste(x, collapse="_"))

dset <- list()
for (group in names(sets)){
  dset1 <- in_silico_perturbation(pset=data_z[, colnames(data_z) == "CMP"], rset=data_z,
                                  up=names(sets[[group]][sets[[group]] > 0]),
                                  dn=names(sets[[group]][sets[[group]] < 0]), strategy="v1", quantile=0.01)
  colnames(dset1) <- rep(paste("GMP", group, sep="_"), ncol(dset1))
  dset2 <- in_silico_perturbation(pset=data_z[, colnames(data_z) == "CMP"], rset=data_z,
                                  up=names(sets[[group]][sets[[group]] < 0]),
                                  dn=names(sets[[group]][sets[[group]] > 0]), strategy="v1", quantile=0.01)
  colnames(dset2) <- rep(paste("MEP", group, sep="_"), ncol(dset2))
  dset[[group]] <- cbind(dset1, dset2)}
dset <- cbind(data_z, do.call(cbind, dset))

gz <- gzfile("../screen/z_screen.csv.gz", "w")
write.table(t(dset), gz, quote=F, row.names=F, col.names=F)
close(gz)
write.table(colnames(dset), file="../screen/annotation_screen.txt", quote=F, row.names=F, col.names=F)

###############################################
##########    prepare screen data    ##########
###############################################
