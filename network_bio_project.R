#BiocManager::install(c("igraph", "RCy3", "gProfileR"))
#BiocManager::install("paxtoolsr")
#BiocManager::install("linkprediction")
library(RCy3)
library(igraph)
library(gProfileR)
library(paxtoolsr)
library(linkprediction)

#cytoscapePing ()
#networkSuid <- getNetworkSuid()

#directory of file
f <- '/home/esther/Documents/Uni/SB/Block7/Network Bio/project/CTL_processed_repoDB.sif'
df <- read.delim(f) # read sif file as tab delimited file into dataframe
gr <- graph_from_data_frame(df, directed = FALSE) # Create graph from dataframe
comps <- components(gr)
largest_comp <- which(comps$membership == which.max(comps$csize))
g <- induced_subgraph(gr, largest_comp)
#plot(g)

iterations <- 50
metrics <- matrix(0, iterations, 4)

for(i in 1:iterations){
  # Returns the type of all edges
  types <- edge_attr(g, "Type")
  #select Drug_Disease interactions
  DD_edges <- which(types=="DD")
  #sample 90% of the DD edges
  
  #plot for validation
  #plot(subg)
  usable_graph <- FALSE
  while(!usable_graph){
    sampled <- sample(DD_edges, 6600, replace = FALSE)
    
    #get set of removed edges
    removed_edge <- setdiff(DD_edges, sampled)
    
    other_edges <- which(types!="DD")
    #combine sampled edges from DD interactions with all others
    combined_edges <- c(other_edges, sampled)
    
    #make subgraph using these edges
    subg <- subgraph.edges(g, combined_edges, delete.vertices = FALSE)
    subg <- as.undirected(subg)
    comps <- components(subg)
    message("After sampling the graph has ",comps$no, " components")
    if(comps$no > 1){
      message("Resampling!")
    }
    else{
      usable_graph <- TRUE
    }
  }
  
  #Predicting new links using common neighbors method
  res_edgelist <- proxfun(subg, method = "katz", value = "edgelist")
  
  
  #Analyse results
  #We want to get the new drug - indication links
  new_edges <- res_edgelist[which(res_edgelist$value >0),]
  select_types_from <- vertex_attr(subg, "name", new_edges$from)
  select_types_to <- vertex_attr(subg, "name", new_edges$to)
  
  interactions <- cbind(select_types_from, select_types_to)
  index_from_drugs <- grep("DB+", select_types_from)
  index_to_indication <- grep("C+", select_types_to)
  drug_indication_index <- intersect(index_from_drugs, index_to_indication) # indices where both these conditions hold
  
  index_from_indication <- grep("C+", select_types_from)
  index_to_drugs <- grep("DB+", select_types_to)
  indication_drug_index <- intersect(index_from_indication, index_to_drugs)
  
  length(indication_drug_index)
  length(drug_indication_index)
  
  drug_indication_links <- new_edges[c(indication_drug_index, drug_indication_index),]
  
  # Remove the edges that were already part of the subnetwork (Why would proxfun return these???)
  from_to <- cbind(drug_indication_links$from, drug_indication_links$to)
  vector_vertices <- as.vector(t(from_to))
  found_edges <- get.edge.ids(subg, vector_vertices, error = FALSE)
  # Filter drug_indication_links to remove all the edges that are already in the subnetwork
  drug_indication_links <- drug_indication_links[found_edges == 0, ]
  
  # select predicted edges with a threshold
  threshold <-2.8e-7
  # links that are above the threshold (new edges)
  drug_indication_links_sign <- drug_indication_links[drug_indication_links$value >= threshold,]
  # links that are below the threshold (deemed not to be edges)
  negatives <- drug_indication_links[drug_indication_links$value < threshold,]
  
  # Remove duplicates. Network is undirected so we only want the links in one direction
  # We do this after thresholding because the score can be different from drug to indication vs. indication to drug
  is_duplicate <- !duplicated(cbind(pmin(drug_indication_links_sign$from, drug_indication_links_sign$to), pmax(drug_indication_links_sign$from, drug_indication_links_sign$to)))
  drug_indication_links_sign <- drug_indication_links_sign[is_duplicate,]
  
  is_duplicate <- !duplicated(cbind(pmin(negatives$from, negatives$to), pmax(negatives$from, negatives$to)))
  negatives <- negatives[is_duplicate,]
  
  #validate links we predicted
  from_to <- cbind(drug_indication_links_sign$from, drug_indication_links_sign$to)
  vector_vertices <- as.vector(t(from_to))
  # check if we can find the edges we predicted
  # We can use found vertex ids in the original network, because we did not remove any vertices (only edges)
  
  found_edges <- get.edge.ids(g, vector_vertices, error = FALSE)
  
  #validate the links we did NOT predict
  negative_from_to <- cbind(negatives$from, negatives$to)
  negative_vector_vertices <- as.vector(t(negative_from_to))
  found_negative_edges <- get.edge.ids(g, negative_vector_vertices, error = FALSE)
  
  # calculate metrics
  # True Positives (Predicted edges that were part of the original network)
  tp <- sum(found_edges != 0)
  # False Positives (Predicted edges that were not part of the original network)
  fp <- sum(found_edges == 0)
  # True negatives (Edges that were not part of the original network and not predicted as such)
  tn <- sum(found_negative_edges == 0) # should be equal to 0, because they should not exist in g
  # False negatives (edges that were part of the original network, but not predicted as edges)
  fn <- sum(found_negative_edges != 0) # should not be 0, because these exist in g
  if(tp != sum(removed_edge %in% found_edges)){
    message("Something wrong!")
  }
  
  message("tp=", tp, " fp=", fp, " tn=", tn, " fn=", fn)
  metrics[i, ] <- c(tp, fp, tn, fn)
  message("Finished iteration ", i, "/", iterations)
  i <- i + 1
}

metrics <- as.data.frame(metrics)
colnames(metrics) <- c("TP", "FP", "TN", "FN")



#vertex_attr(g, )





######

#setwd('/home/esther/Documents/Uni/SB/Block7/Network Bio/project')
#g <- loadSifInIgraph('CTL_processed_repoDB.sif')

# lists all attribute names of vertices
#vertex_attr_names(g)
# Lists all attribute names of edges
#edge_attr_names(g)
#vertices_sampled <- vertex_attr(sampled)

# Returns the type of only the 3 first edges
#select_types <- edge_attr(g, "Type", c(1,2,3))

#Add found edges to graph and validate them
#DI_links_idx <- as.vector(t(drug_indication_links))
#ID_links_idx <- as.vector(t(indication_drugs_links))

#res_graph <- add_edges(subg, new_edges_vector, attr = 'new')
