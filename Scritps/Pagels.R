#Pagels lambda assignment
#Simulate a tree
library(ape)
library(phytools)
set.seed(132)
#Generate a random tree
tree <- rtree(n = 4)
#Generates the random trait values
traits <- rnorm(4)
#Plot tree and visualize inner nodes
plot(tree)
nodelabels()
#Determine branch length for calculating matrix
branch_lengths <- tree$edge.length
print(branch_lengths)
# Print the edge matrix so you can figure out where those branch lengths are
edge_matrix <- tree$edge
print(edge_matrix)
#Write matrix by hand 

#Get pagels lambda calculation is correct using phylosig
lambda_result <- phylosig(tree, traits, method = "lambda")
print(lambda_result)
