##########################################################
####   IDS564 Project: Analysis of FLixster Network   ####
####                     Rahul Shah                   ####
##########################################################

## Setting Working Directory 
dir_path <- getwd()
setwd(dir_path)

## Clearing everything out of memory
rm(list=ls()) 

## Load package
library(igraph)

## Reading edge and node data
el = read.csv("edges.csv", header=F)
nl = read.csv("nodes.csv", header=F)
colnames(el) <- c("Source","Target")
colnames(nl) <- c("Nodes")
# write.csv(el, "edges_gephi.csv", row.names=FALSE)
# write.csv(nl, "nodes_gephi.csv", row.names=FALSE)
nl1 = nl

class(el)
# ---
# [1] "data.frame"
# ---
# Describe the data frame

str(el)

## Creating the undirected graph object for Flixster network
g_Flixster = graph.data.frame(el, directed = FALSE, vertices= nl)

## Edges
ecount(g_Flixster)
# ---
# [1] 9197337
# ---

## Vertices
vcount(g_Flixster)
# ---
# [1] 2523386
# ---

## Generating an induced subgraph of the Flixster network
nl_subset = data.frame(V1 = nl[1:5000,])
g_Flixter_subgraph = induced_subgraph(g_Flixster, nl_subset[["V1"]])

## Edges
ecount(g_Flixter_subgraph)
# ---
# [1] 49009
# ---

## Vertices
vcount(g_Flixter_subgraph)
# ---
# [1] 5000
# ---

## Exporting the subgraph for use in Gephi
el_subset = data.frame(get.edgelist(g_Flixter_subgraph))
colnames(el_subset) <- c("Source","Target")
colnames(nl_subset) <- c("Nodes")
write.csv(el_subset, "edges_subgraph.csv", row.names=FALSE)
write.csv(nl_subset, "nodes_subgraph.csv", row.names=FALSE)

plot(g_Flixter_subgraph, layout=layout.fruchterman.reingold(g_Flixter_subgraph))
# write.csv(edge_frame_subset,file = "Subset_Edges.csv")

## Check whether Self_loops exist, as do multiple edges
is.simple(g_Flixter_subgraph)
#Is it a simple graph? No!
# ---
#[1] FALSE
# ---

## Simplifying the network graph 
g_Flixster_simpl = simplify(g_Flixter_subgraph) 
is.simple(g_Flixster_simpl)
#Is it a simple graph? Yes!
# ---
#[1] TRUE
# ---

## Summarize the graph structure
summary(g_Flixster_simpl)

## Edges
ecount(g_Flixster_simpl)
# ---
# [1] 24507
# ---

## Vertices
vcount(g_Flixster_simpl)
# ---
# [1] 5000
# ---

## Exporting the simplified subgraph for use in Gephi
el_subset_simpl = data.frame(get.edgelist(g_Flixster_simpl))
colnames(el_subset_simpl) <- c("Source","Target")
write.csv(el_subset_simpl, "edges_subgraph_simpl.csv", row.names=FALSE)


############################
#### Graph connectivity ####
############################
# We can look at the connectivity of the graph 
# by seeing the neighbors of some selected nodes
neighbors(g_Flixster_simpl, v=c('500'))
# ---
## + 11/5000 vertices, named, from ccab81e:
##  [1] 5    20   316  407  591  1543 1563 3812 4039 4139 4319
# ---

neighbors(g_Flixster_simpl, v=c('4925'))
# ---
## + 2/5000 vertices, named, from ccab81e:
## [1] 37   2210
# ---


######################
#### Node Degrees ####
######################
# All Degrees
all.deg <- degree(g_Flixster_simpl, v=V(g_Flixster_simpl), mode="all")
table(all.deg)

# Overall Max Degree Node
max(all.deg)
# ---
# [1] 731
# ---

V(g_Flixster_simpl)$name[degree(g_Flixster_simpl)==max(degree(g_Flixster_simpl))] # gives the corresponding user
# ---
# [1] "21"
# ---

# Mean Degree
# table(degree(g_Flixster_simpl, mode="all"))
mean(degree(g_Flixster_simpl), mode="all")
# ---
# [1] 9.8028
# ---


##########################################
# Vertex Strength Distribution Histogram #
##########################################
# par(mfrow=c(2,1))
hist(graph.strength(g_Flixster_simpl, mode='all'), col="lightblue",
     xlab="Vertex Strength", ylab="Frequency", 
     main="Vertex Strength Distribution",
     breaks = 100)


###############################
# Log-log degree distribution #
###############################
# par(mfrow=c(1, 2))
d.net <- degree(g_Flixster_simpl, mode='all')
dd.net <- degree.distribution(g_Flixster_simpl, mode='all')
d <- 1:max(d.net)-1
ind <- (dd.net != 0)
plot(d[ind], dd.net[ind], log="xy", col="blue",
     xlab=c("Log-Degree"), ylab=c("Log-Intensity"),
     main="Log-Log Degree Distribution")


###############################################
# Average neighbor degree versus vertex degree#
###############################################
par(mfrow=c(1, 1))
d.net <- degree(g_Flixster_simpl, mode='all')
a.nn.deg <- graph.knn(g_Flixster_simpl,V(g_Flixster_simpl))$knn
plot(d.net, a.nn.deg, log="xy", 
     col="goldenrod", xlab=c("Log Vertex Degree"),
     ylab=c("Log Average Neighbor Degree"))

assortativity(g_Flixster_simpl, types1 = V(g_Flixster_simpl), directed = TRUE)
assortativity_nominal(g_Flixster_simpl, types=V(g_Flixster_simpl), directed = FALSE)
assortativity_degree(g_Flixster_simpl, directed = FALSE)

####################################
# Connectivity Measures of Network #
#################################### 
## Network Density
graph.density(g_Flixster_simpl)
# ---
## [1] 0.001960952
# ---

## Average Path Length
mean_distance(g_Flixster_simpl) #average.path.length(g_Flixster_simpl, directed=TRUE)
# ---
## [1] 3.488269
# ---

## Number of weakly connected components
is.connected(g_Flixster_simpl, mode="weak")
# ---
## [1] TRUE
# ---
g_Flixster_simpl.wcc <- clusters(g_Flixster_simpl, mode="weak")
table(g_Flixster_simpl.wcc$csize)
# ---
## 5000 
##    1
# ---

## Number of strongly connected components
is.connected(g_Flixster_simpl, mode="strong")
# ---
## [1] TRUE
# ---
g_Flixster_simpl.scc <- clusters(g_Flixster_simpl, mode="strong")
table(g_Flixster_simpl.scc$csize)
# ---
## 5000 
##    1
# ---

# ########################
# # Reciprocal Relations #
# ########################
# reciprocity(g_Flixster_simpl)
# dyad_census(g_Flixster_simpl)
# which_mutual(g_Flixster_simpl)
# sum(which_mutual(g_Flixster_simpl))/2 == dyad_census(g_Flixster_simpl)$mut
# # Table of vertexes in mutual edges
# table(unlist(strsplit(attr(E(g_Flixster_simpl)[which_mutual(g_Flixster_simpl)],"vnames"),"\\|")))/2

####################
# Network Diameter #
####################
diameter(g_Flixster_simpl)
# ---
## [1] 4
# ---

##############
# Clustering #
##############
# Global clustering coefficient
transitivity(g_Flixster_simpl)
# ---
## [1] 0.07100749
# --- 

#################
# Clique Census #
#################
# Clique structure: 5 cliques of size 5, 39 cliques of size 4, 335 triangles
clique.number(g_Flixster_simpl)
# ---
## [1] 16
# --- 

table(sapply(cliques(g_Flixster_simpl), length))
# ---
##    1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16 
## 5000  24507  31181  51980  88648 138202 184421 201254 175724 121081  64833  26381   7875   1624    206     12
# ---

table(sapply(maximal.cliques(g_Flixster_simpl), length))
# ---
##    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
## 8275 4478 1925 1405  889  464  250  125  103  130   53   77   34   30   12
# ---

# cliques(g_Flixster_simpl)[sapply(cliques(g_Flixster_simpl), length) == 16]

# A <- get.adjacency(g_Flixster_simpl, sparse=FALSE)


################################
# Local Measures: Centralities #
################################

# Embeddedness/ inverse of structural hole access (see Burt 2004)
constraints_Flixster <- round(constraint(g_Flixster_simpl, nodes=V(g_Flixster_simpl)), digits=4)
# Degree centrality
degree_Flixster <- degree(g_Flixster_simpl)
# Node betweenness
nodebetweens_Flixster <- round(betweenness(g_Flixster_simpl, v=V(g_Flixster_simpl), directed = TRUE, nobigint =TRUE, normalized = FALSE))
# Edge betwenness
edgebetweens_Flixster<-edge.betweenness(g_Flixster_simpl, e=E(g_Flixster_simpl), directed = TRUE)
# Local clustering coefficients
clustering_Flixster <- transitivity(g_Flixster_simpl, type="local", vids=V(g_Flixster_simpl)) 

# # Plots 1 and 2: Can run them together
# par(mfrow=c(1, 1))
# edge_frame<-data.frame(edgebetweens_Flixster, num_weight, inv_weight)
# a_edge<-aggregate(edgebetweens_Flixster ~ inv_weight, data=edge_frame, mean)
# plot(a_edge, col="blue", log="xy", xlab="Weight of edge", ylab="Average Betweenness of edges")
# node_frame<-data.frame(nodebetweens_Flixster, constraints_Flixster, clustering_Flixster, degree_Flixster)
# a_node<-aggregate(nodebetweens_Flixster ~ clustering_Flixster, data=node_frame, mean)
# plot(a_node, col="blue", log="xy", xlab="Clustering", ylab="Average Betweenness of nodes")

# Plot set 2: Four plots 
par(mfrow=c(2, 2))

node_frame<-data.frame(nodebetweens_Flixster, constraints_Flixster, clustering_Flixster, degree_Flixster)
a_node<-aggregate(nodebetweens_Flixster ~ clustering_Flixster, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Clustering", ylab="Average Betweenness of nodes")

a_node<-aggregate(nodebetweens_Flixster ~ degree_Flixster, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Betweenness of nodes")

# a_edge<-aggregate(edgebetweens_Flixster ~ num_weight, data=edge_frame, mean)
# plot(a_edge, col="blue", log="xy", xlab="Weight of edge", ylab="Average Betweenness of edges")

a_node<-aggregate(clustering_Flixster ~ degree_Flixster, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Clustering")

a_node<-aggregate(constraints_Flixster ~ degree_Flixster, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Constraint (Embeddedness)")

# Closeness Centrality
close_Flixster <- closeness(g_Flixster_simpl)
# Eigen Centrality
eig_Flixster <- evcent(g_Flixster_simpl)$vector

# # Hub and Authority Scores
# hub_Flixster <- hub.score(g_Flixster_simpl, weights=inv_weight)$vector
# auth_Flixster <- authority.score(g_Flixster_simpl, weights=inv_weight)$vector
# head(sort(hub_Flixster, decreasing=TRUE))
# head(sort(auth_Flixster, decreasing=TRUE))

centralities <- cbind(degree_Flixster, nodebetweens_Flixster, edgebetweens_Flixster, close_Flixster, eig_Flixster)
cor.matrix = round(cor(centralities), 4)
write.csv(cor.matrix, "cor_matrix.csv")


#############################
###  Community Detection  ###
#############################
# "The following code snippet performs a Wilcoxon rank-sum test on the "internal" and "external"
# degrees of a community in order to quantify its significance. Let us call the edges within a 
# community "internal" and the edges connecting the vertices of a community with the rest of the graph "external".
# The null hypothesis of the test is that there is no difference between the number of "internal" and "external" edges 
# incident to a vertex of the community. More internal than external edges show that the community is significant; less 
# internal than external edges show that the community is in fact an "anti-community". The p-value of the test performed by 
# this function will be close to zero in both cases; the value of the test statistic tells us whether we have a community or an anti-community."

community.significance.test <- function(graph, vs, ...) {
  if (is.directed(graph)) stop("This method requires an undirected graph")
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  # Total degree among nodes in the vs list, minus the degree within the subgraph 
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}

#####################################################################
####### Community detection using the Fast Greedy Algorithm #########
#####################################################################

flixster_comm_fast <- fastgreedy.community(g_Flixster_simpl)
c.m.fast <- membership(flixster_comm_fast)
table(c.m.fast, useNA = c("no"))

# Testing community significance for various communities
v_comp1 <- V(g_Flixster_simpl)[c.m.fast==1]
v_comp2 <- V(g_Flixster_simpl)[c.m.fast==2]
v_comp3 <- V(g_Flixster_simpl)[c.m.fast==3]
v_comp4 <- V(g_Flixster_simpl)[c.m.fast==4]
v_comp5 <- V(g_Flixster_simpl)[c.m.fast==5]
v_comp6 <- V(g_Flixster_simpl)[c.m.fast==6]
v_comp7 <- V(g_Flixster_simpl)[c.m.fast==7]
v_comp8 <- V(g_Flixster_simpl)[c.m.fast==8]
v_comp9 <- V(g_Flixster_simpl)[c.m.fast==9]
v_comp10 <- V(g_Flixster_simpl)[c.m.fast==10]
v_comp11 <- V(g_Flixster_simpl)[c.m.fast==11]
v_comp12 <- V(g_Flixster_simpl)[c.m.fast==12]
v_comp13 <- V(g_Flixster_simpl)[c.m.fast==13]
v_comp14 <- V(g_Flixster_simpl)[c.m.fast==14]
v_comp15 <- V(g_Flixster_simpl)[c.m.fast==15]
v_comp16 <- V(g_Flixster_simpl)[c.m.fast==16]
v_comp17 <- V(g_Flixster_simpl)[c.m.fast==17]
v_comp18 <- V(g_Flixster_simpl)[c.m.fast==18]
community.significance.test(g_Flixster_simpl, v_comp1) # yes: W = 1081100, p-value < 2.2e-16
community.significance.test(g_Flixster_simpl, v_comp2) # yes: W = 57082, p-value = 0.0007916
community.significance.test(g_Flixster_simpl, v_comp3) # yes: W = 17565, p-value = 2.922e-08
community.significance.test(g_Flixster_simpl, v_comp4) # yes: W = 1247900, p-value < 2.2e-16
community.significance.test(g_Flixster_simpl, v_comp5) # yes: W = 27517, p-value < 2.2e-16
community.significance.test(g_Flixster_simpl, v_comp6) # yes: W = 285970, p-value < 2.2e-16
community.significance.test(g_Flixster_simpl, v_comp7) # yes: W = 171300, p-value < 2.2e-16
community.significance.test(g_Flixster_simpl, v_comp8) # yes: W = 23764, p-value < 2.2e-16
community.significance.test(g_Flixster_simpl, v_comp9) # yes: W = 37986, p-value < 2.2e-16
community.significance.test(g_Flixster_simpl, v_comp10) # yes: W = 45776, p-value = 2.04e-13
community.significance.test(g_Flixster_simpl, v_comp11) # yes: W = 2903.5, p-value = 8.926e-12
community.significance.test(g_Flixster_simpl, v_comp12) # yes: W = 1177, p-value = 4.212e-08
community.significance.test(g_Flixster_simpl, v_comp13) # yes: W = 1573.5, p-value = 4.606e-09
community.significance.test(g_Flixster_simpl, v_comp14) # no: W = 39.5, p-value = 0.9627
community.significance.test(g_Flixster_simpl, v_comp15) # yes: W = 59.5, p-value = 0.002333
community.significance.test(g_Flixster_simpl, v_comp16) # yes: W = 100.5, p-value = 0.004454
community.significance.test(g_Flixster_simpl, v_comp17) # no: W = 27, p-value = 0.07401
community.significance.test(g_Flixster_simpl, v_comp18) # yes: W = 87, p-value = 0.00348

# -------------------------------------------------------------------------
## c.m  | Test Stat, W | p-value  | Significant | Community/Anti-Community
## -----|--------------|----------|-------------|-------------------------- 
##   1  |   1512.0     | 0.002056 |     Yes     |       Community
##   2  |    321.5     | 0.494    |      No     |         ----
##   3  |   1321.0     | 0.003201 |     Yes     |       Community
##   4  |    276.0     | 0.8079   |      No     |         ----
##   5  |   1911.5     | 4.356e-07|     Yes     |       Community
##   6  |    279.5     | 0.5248   |      No     |         ----
##   7  |    273.5     | 0.2351   |      No     |         ----
# -------------------------------------------------------------------------

par(mfrow=c(1, 1))
plot(flixster_comm_fast, g_Flixster_simpl, vertex.label= NA, vertex.size=4)
# plot(flixster_comm_fast, g_Flixster_simpl, layout=layout.lgl, vertex.label= NA, vertex.size = 4)
# plot(flixster_comm_fast, g_Flixster_simpl, layout = layout.fruchterman.reingold(g_Flixster_simpl), vertex.label= NA, vertex.size=4)
# plot(flixster_comm_fast, g_Flixster_simpl, layout = layout.kamada.kawai(g_Flixster_simpl), vertex.label= NA, vertex.size=4)
# # plot(flixster_comm_fast, g_Flixster_simpl, layout = layout.forceatlas2(g_Flixster_simpl), vertex.label= NA, vertex.size=4)
# 
# install.packages("devtools")
# if (!require("ForceAtlas2")) devtools::install_github("analyxcompany/ForceAtlas2")
# library("ForceAtlas2")
# layout.forceatlas2(g_Flixster_simpl, iterations=10000, plotstep=500)

sub_comp1 <- induced.subgraph(g_Flixster_simpl, v=v_comp4) # 1252
sub_comp2 <- induced.subgraph(g_Flixster_simpl, v=v_comp7) # 500
sub_comp3 <- induced.subgraph(g_Flixster_simpl, v=v_comp8) # 178
sub_comp4 <- induced.subgraph(g_Flixster_simpl, v=v_comp11) # 58
sub_comp5 <- induced.subgraph(g_Flixster_simpl, v=v_comp12) # 37
sub_comp6 <- induced.subgraph(g_Flixster_simpl, v=v_comp16) # 11

plot(sub_comp1, layout=layout.fruchterman.reingold(sub_comp1), vertex.label= NA, vertex.size = 6)
plot(sub_comp2, layout=layout.kamada.kawai(sub_comp2), vertex.label= NA, vertex.size = 6)
plot(sub_comp3, layout=layout.kamada.kawai(sub_comp3), vertex.label= NA, vertex.size = 6)
plot(sub_comp4, layout=layout.kamada.kawai(sub_comp4), vertex.label= NA, vertex.size = 6)
plot(sub_comp5, layout=layout.kamada.kawai(sub_comp5), vertex.label= NA, vertex.size = 6)
plot(sub_comp6, layout=layout.kamada.kawai(sub_comp6), vertex.label= NA, vertex.size = 6)

##################################################################
####### Community detection using the Walktrap Algorithm #########
##################################################################

flixster_comm_walktrap <- walktrap.community(g_Flixster_simpl)
c.m.walk <- membership(flixster_comm_walktrap)
table(c.m.walk, useNA = c("no"))

c.m.walk.table = table(c.m.walk, useNA = c("no"))
length(c.m.walk.table[c.m.walk.table>20])

# Testing Community Significance
v_comp1_w <- V(g_Flixster_simpl)[c.m.walk==1]
v_comp2_w <- V(g_Flixster_simpl)[c.m.walk==2]
v_comp3_w <- V(g_Flixster_simpl)[c.m.walk==3]
v_comp4_w <- V(g_Flixster_simpl)[c.m.walk==4]
v_comp5_w <- V(g_Flixster_simpl)[c.m.walk==5]
v_comp6_w <- V(g_Flixster_simpl)[c.m.walk==6]
v_comp7_w <- V(g_Flixster_simpl)[c.m.walk==7]
v_comp8_w <- V(g_Flixster_simpl)[c.m.walk==8]
v_comp9_w <- V(g_Flixster_simpl)[c.m.walk==9]

community.significance.test(g_Flixster_simpl, v_comp1_w)
community.significance.test(g_Flixster_simpl, v_comp2_w)
community.significance.test(g_Flixster_simpl, v_comp3_w)
community.significance.test(g_Flixster_simpl, v_comp4_w)
community.significance.test(g_Flixster_simpl, v_comp5_w)
community.significance.test(g_Flixster_simpl, v_comp6_w)
community.significance.test(g_Flixster_simpl, v_comp7_w)
community.significance.test(g_Flixster_simpl, v_comp8_w)
community.significance.test(g_Flixster_simpl, v_comp9_w)

# Network Graph Plot Showing Community Structure Using Walktrap Algorithm
plot(flixster_comm_walktrap, g_Flixster_simpl, vertex.label= NA, vertex.size=4)
# plot(flixster_comm_walktrap, g_Flixster_simpl, layout = layout.fruchterman.reingold(g_Flixster_simpl), vertex.label= NA, vertex.size=4)
# plot(flixster_comm_walktrap, g_Flixster_simpl, layout = layout.kamada.kawai(g_Flixster_simpl), vertex.label= NA, vertex.size=4)

#################################################################
####### Community detection using the Spinglass Algorithm #######
#################################################################

set.seed(1234)
flixster_comm_spinglass <- spinglass.community(g_Flixster_simpl)
c.m.spinglass <- membership(flixster_comm_spinglass)
table(c.m.spinglass, useNA = c("no"))

# Testing Community Significance
v_comp1_s <- V(g_Flixster_simpl)[c.m.spinglass==1]
v_comp2_s <- V(g_Flixster_simpl)[c.m.spinglass==2]
v_comp3_s <- V(g_Flixster_simpl)[c.m.spinglass==3]
v_comp4_s <- V(g_Flixster_simpl)[c.m.spinglass==4]
v_comp5_s <- V(g_Flixster_simpl)[c.m.spinglass==5]
v_comp6_s <- V(g_Flixster_simpl)[c.m.spinglass==6]
v_comp7_s <- V(g_Flixster_simpl)[c.m.spinglass==7]

community.significance.test(g_Flixster_simpl, v_comp1_s)
community.significance.test(g_Flixster_simpl, v_comp2_s)
community.significance.test(g_Flixster_simpl, v_comp3_s)
community.significance.test(g_Flixster_simpl, v_comp4_s)
community.significance.test(g_Flixster_simpl, v_comp5_s)
community.significance.test(g_Flixster_simpl, v_comp6_s)
community.significance.test(g_Flixster_simpl, v_comp7_s)

# Network Graph Plot Showing Community Structure Using Spinglass Algorithm
plot(flixster_comm_spinglass, g_Flixster_simpl, vertex.label= NA, vertex.size=4)
# plot(flixster_comm_spinglass, g_Flixster_simpl, layout = layout.fruchterman.reingold(g_Flixster_simpl), vertex.label= NA, vertex.size=4)
# plot(flixster_comm_spinglass, g_Flixster_simpl, layout = layout.kamada.kawai(g_Flixster_simpl), vertex.label= NA, vertex.size=4)

###########################################################################
####### Community detection using the Label Propagation Algorithm #########
###########################################################################

set.seed(2)
flixster_comm_label <- label.propagation.community(g_Flixster_simpl)
c.m.label <- membership(flixster_comm_label)
table(c.m.label, useNA = c("no"))

# Testing Community Significance
v_comp1_lp <- V(g_Flixster_simpl)[c.m.label==1]
v_comp2_lp <- V(g_Flixster_simpl)[c.m.label==2]
v_comp3_lp <- V(g_Flixster_simpl)[c.m.label==3]
v_comp4_lp <- V(g_Flixster_simpl)[c.m.label==4]
v_comp5_lp <- V(g_Flixster_simpl)[c.m.label==5]
v_comp6_lp <- V(g_Flixster_simpl)[c.m.label==6]
v_comp7_lp <- V(g_Flixster_simpl)[c.m.label==7]
v_comp8_lp <- V(g_Flixster_simpl)[c.m.label==8]
v_comp9_lp <- V(g_Flixster_simpl)[c.m.label==9]
v_comp10_lp <- V(g_Flixster_simpl)[c.m.label==10]
v_comp11_lp <- V(g_Flixster_simpl)[c.m.label==11]

community.significance.test(g_Flixster_simpl, v_comp1_lp)
community.significance.test(g_Flixster_simpl, v_comp2_lp)
community.significance.test(g_Flixster_simpl, v_comp3_lp)
community.significance.test(g_Flixster_simpl, v_comp4_lp)
community.significance.test(g_Flixster_simpl, v_comp5_lp)
community.significance.test(g_Flixster_simpl, v_comp6_lp)
community.significance.test(g_Flixster_simpl, v_comp7_lp)
community.significance.test(g_Flixster_simpl, v_comp8_lp)
community.significance.test(g_Flixster_simpl, v_comp9_lp)
community.significance.test(g_Flixster_simpl, v_comp10_lp)
community.significance.test(g_Flixster_simpl, v_comp11_lp)

# Network Graph Plot Showing Community Structure Using Label Propagation Algorithm
plot(flixster_comm_label, g_Flixster_simpl, vertex.label= NA, vertex.size=4)
# plot(flixster_comm_label, g_Flixster_simpl, layout = layout.fruchterman.reingold(g_Flixster_simpl), vertex.label= NA, vertex.size=4)
# plot(flixster_comm_label, g_Flixster_simpl, layout = layout.kamada.kawai(g_Flixster_simpl), vertex.label= NA, vertex.size=4)

#####################################################################
####### Community detection using the Girvan-Newman Algorithm #######
#####################################################################
flixster_comm_Girvan_Newman <- edge.betweenness.community(g_Flixster_simpl)
c.m.gm <- membership(flixster_comm_Girvan_Newman)
table(c.m.gm, useNA = c("no"))
length(flixster_comm_Girvan_Newman)
sizes(flixster_comm_Girvan_Newman)

# Network Graph Plot Showing Community Structure Using Girvan-Newman Algorithm
plot(flixster_comm_Girvan_Newman, g_Flixster_simpl, vertex.label= NA, vertex.size=4)
# plot(flixster_comm_Girvan_Newman, g_Flixster_simpl, layout = layout.fruchterman.reingold(g_Flixster_simpl), vertex.label= NA, vertex.size=4)
# plot(flixster_comm_Girvan_Newman, g_Flixster_simpl, layout = layout.kamada.kawai(g_Flixster_simpl), vertex.label= NA, vertex.size=4)


##############################################
# Number of Common Neighbors and Viola Plots #
##############################################
nv <- vcount(g_Flixster_simpl)
ncn <- numeric()
A <- get.adjacency(g_Flixster_simpl)

# plot(g_Flixster_simpl, vertex.size=3, vertex.label=NA)

# Find the number of common neighbors for each pair of nodes in the fblog network
for(i in (1:(nv-1))){
  ni <- neighborhood(g_Flixster_simpl, 1, i)
  nj <- neighborhood(g_Flixster_simpl, 1, (i+1):nv)
  nbhd.ij <- mapply(intersect, ni, nj, SIMPLIFY=FALSE)
  temp <- unlist(lapply(nbhd.ij, length)) - 
    2*A[i, (i+1):nv]
  ncn <- c(ncn, temp)
}

library(vioplot)
Avec <- A[lower.tri(A)]
vioplot(ncn[Avec==0], ncn[Avec==1], 
        names=c("No Edge", "Edge"), col='purple')
title(ylab="Number of Common Neighbors")

# library(ROCR)
# pred <- prediction(ncn, Avec)
# perf <- performance(pred, "auc")
# slot(perf, "y.values")