\name{assign_to_tree}
\alias{assign_to_tree}
\title{Assigns mutations to a tree topology and sets the branch lengths based on the hard assignment of mutations to branches.}

\usage{assign_to_tree(tree, mtr, depth, error_rate = rep(0.01, dim(mtr)[2]), 
    maxits = 5)}
\arguments{
  \item{tree}{Tree. Should already include the zeros outgroup..}
  \item{mtr}{Mutant read count matrix. Each row is a mutation and each colomn is a clonal sample}
  \item{depth}{Depth count matrix. Each row is a mutation and each colomn is a clonal sample}
  \item{error_rate}{Error rate. This is intended as a coarse error measure and should not take more than 4 distinct values. Useful if a sample is mildly contaminated error->0.1.}
  \item{maxits}{
}
}

\value{Returns a list containing:
tree: the adjusted tree
and a data.frame "summary" that is aligned with the input matrices and maps the mutations to branches:
summary$edge_ml : Specifies the index of the branch in the tree edge matrix.  
summary$pval :  A heuristic pvalue assessing the hypothesis that the mutation is consistent with the provided tree topology.   }

\author{Nick Williams}





