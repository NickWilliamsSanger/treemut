### R code from vignette source 'rclonetree.Rnw'

###################################################
### code chunk number 1: rclonetree.Rnw:15-20
###################################################
  options(continue="  ")
options(digits=4)
set.seed(323123)
require("rclonetree")
require("xtable")


###################################################
### code chunk number 2: rclonetree.Rnw:28-29
###################################################
 tree=generate_random_tree(50)


###################################################
### code chunk number 3: splot
###################################################
plot(tree)


###################################################
### code chunk number 4: rclonetree.Rnw:39-40
###################################################
plot(tree)


###################################################
### code chunk number 5: rclonetree.Rnw:46-47 (eval = FALSE)
###################################################
##  tree=bind.tree(tree,read.tree(text="(zeros:0);"))


###################################################
### code chunk number 6: rclonetree.Rnw:51-54
###################################################
  df=reconstruct_genotype_summary(tree)
cat(df$samples,"\n")
head(df$df[, 1:3])


###################################################
### code chunk number 7: rclonetree.Rnw:68-72
###################################################
 simdat=simulate_reads_from_tree(df,12)
head(simdat$mtr)
head(simdat$depth)
print(simdat$p.error)


###################################################
### code chunk number 8: rclonetree.Rnw:76-77
###################################################
 res=assign_to_tree(tree,simdat$mtr,simdat$depth,error_rate=simdat$p.error)


###################################################
### code chunk number 9: splot2
###################################################
tree_estimated=res$tree
par(mfcol=c(1,2))
plot(ladderize(tree,right=TRUE),cex=0.5)
plot(ladderize(tree_estimated,right=TRUE),cex=0.5)


###################################################
### code chunk number 10: rclonetree.Rnw:89-90
###################################################
tree_estimated=res$tree
par(mfcol=c(1,2))
plot(ladderize(tree,right=TRUE),cex=0.5)
plot(ladderize(tree_estimated,right=TRUE),cex=0.5)


###################################################
### code chunk number 11: splot3
###################################################
sim=list(edge_length_orig=df$df$edge_length,
         edge_length_inferred=res$df$df$edge_length,
         expected_edge_length_inferred=res$df$df$expected_edge_length,
         edge_idx_orig=simdat$edge,
         edge_idx_ml=res$summary$edge_ml)
plot_sim_result(sim)


###################################################
### code chunk number 12: rclonetree.Rnw:106-107
###################################################
sim=list(edge_length_orig=df$df$edge_length,
         edge_length_inferred=res$df$df$edge_length,
         expected_edge_length_inferred=res$df$df$expected_edge_length,
         edge_idx_orig=simdat$edge,
         edge_idx_ml=res$summary$edge_ml)
plot_sim_result(sim)


