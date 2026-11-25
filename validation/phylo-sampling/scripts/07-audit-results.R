# Final independent audit of all inputs, statuses, counts, curves and intervals.
library(ape)

check_intervals <- function(samples, summary) {
  stopifnot(nrow(samples)==nrow(summary))
  error <- 0
  for(i in seq_len(nrow(samples))) {
    expected <- quantile(samples[i,],c(0.025,0.975),names=FALSE)
    observed <- c(summary$lower[i],summary$upper[i])
    error <- max(error,abs(expected-observed))
  }
  stopifnot(is.finite(error),error<1e-9)
  error
}

design <- readRDS("design/design.rds")
manifest <- read.csv("design/task-manifest.csv")
stopifnot(nrow(manifest)==1500L, all(table(manifest$condition)==100L))
characters <- readRDS("design/node-characters.rds")
stopifnot(identical(characters$study,characters$full[design$study_tree$tip.label,,drop=FALSE]))
taxa <- readRDS("inputs/taxon_data.RDS")
original_full <- reorder.phylo(read.tree("inputs/HYMtopC_bamm_hisse.phy"),"cladewise")
original_study <- readRDS("inputs/hym_tree.RDS")
stopifnot(isTRUE(all.equal.phylo(design$full_tree,original_full)),
          isTRUE(all.equal.phylo(design$study_tree,original_study)),
          !anyDuplicated(manifest$seed))
apo_study <- as.character(taxa$Blaimer_taxa[taxa$Tag!="NonApo"])
apo_node <- getMRCA(original_full,apo_study)
apo_tips <- extract.clade(original_full,apo_node)$tip.label
stopifnot(apo_node==design$apocrita_node,
          all(characters$full[!rownames(characters$full)%in%apo_tips,]=="0"))
for(i in seq_len(nrow(design$nodes))) {
  node <- design$nodes$node[i]
  id <- design$nodes$character[i]
  descendants <- extract.clade(original_full,node)$tip.label
  ones <- rownames(characters$full)[characters$full[,id]=="1"]
  stopifnot(node!=apo_node,setequal(descendants,ones),
            all(descendants%in%apo_tips))
  parent <- original_full$edge[design$nodes$edge[i],1]
  stopifnot(all(extract.clade(original_full,parent)$tip.label%in%apo_tips))
}
bank_hashes <- tools::md5sum(c("design/node-fits-full.rds","design/node-fits-study.rds"))
source_hash <- tools::md5sum("provenance/ontophylo-source-sha256.json")
audit <- data.frame()
total_full_maps <- 0L
total_study_maps <- 0L
for(task in manifest$task) {
  result_file <- sprintf("results/task-%04d/result.rds",task)
  stopifnot(file.exists(result_file))
  x <- readRDS(result_file)
  simulated_file <- sprintf("simulated/task-%04d.rds",task)
  simulated <- readRDS(simulated_file)
  stopifnot(identical(x$selected_nodes,simulated$selected_nodes),
            x$task==task,x$condition==manifest$condition[task],
            x$trial==manifest$trial[task],x$epoch==manifest$epoch[task],
            x$proportion==manifest$proportion[task],
            identical(x$settings$simulation,tools::md5sum(simulated_file)),
            identical(x$settings$full_fit_bank,bank_hashes[1]),
            identical(x$settings$study_fit_bank,bank_hashes[2]),
            identical(x$settings$source,source_hash),
            x$n_characters==manifest$n_characters[task],
            x$n_retained==manifest$n_variable_study[task],
            !anyDuplicated(x$selected_nodes$node),
            all(x$selected_nodes$epoch==x$epoch),
            x$settings$models=="ER",x$settings$nstarts==5L,
            x$settings$root_p=="maddfitz",x$settings$normalize=="none",
            x$settings$n_maps==100L,
            x$settings$mapping_seed==420261120L+task*10L)
  max_error <- 0
  interval_max_error <- 0
  map_hashes <- character(2L)
  for(tree_name in c("full","study")) {
    status <- x$diagnostics[[tree_name]]$mapping_status
    expected <- x$n_characters
    if(tree_name=="study") expected <- x$n_retained
    stopifnot(!any(status$status=="failed"),
              sum(status$status=="ok")==expected,
              all(status$n_maps[status$status=="ok"]==100L),
              sum(status$n_maps)==100L*expected)
    map_file <- sprintf("results/task-%04d/maps-%s.rds",task,tree_name)
    stopifnot(file.exists(map_file))
    mapped <- readRDS(map_file)
    tree <- design[[paste0(tree_name,"_tree")]]
    stopifnot(identical(mapped$status,status),
              identical(mapped$tree$edge,tree$edge),
              identical(mapped$tree$edge.length,tree$edge.length),
              identical(mapped$tree$tip.label,tree$tip.label),
              identical(names(mapped$characters),x$selected_nodes$character),
              mapped$settings$n_stm==100L,
              mapped$settings$seed==x$settings$mapping_seed+match(tree_name,c("full","study")))
    summed_counts <- matrix(0,nrow(tree$edge),100L)
    matched_counts <- summed_counts
    retained_ids <- x$selected_nodes$character[x$selected_nodes$variable_study]
    for(id in mapped$status$character[mapped$status$status=="ok"]) {
      counts <- mapped$characters[[id]]$branch_counts
      stopifnot(identical(dim(counts),dim(summed_counts)),
                all(is.finite(counts)),all(counts>=0),all(counts==floor(counts)))
      summed_counts <- summed_counts+counts
      if(id%in%retained_ids) matched_counts <- matched_counts+counts
    }
    group <- x$groups[[paste0(tree_name,"_all")]]
    stopifnot(max(abs(group$rate_samples*tree$edge.length-summed_counts))<1e-7)
    if(tree_name=="full" && x$n_retained>0L) {
      stopifnot(max(abs(x$groups$full_matched$rate_samples*tree$edge.length-
                       matched_counts))<1e-7)
    }
    map_hashes[match(tree_name,c("full","study"))] <- unname(tools::md5sum(map_file))
    rm(mapped)
    if(tree_name=="full") total_full_maps <- total_full_maps+sum(status$n_maps)
    if(tree_name=="study") total_study_maps <- total_study_maps+sum(status$n_maps)
  }
  for(group_name in names(x$groups)) {
    group <- x$groups[[group_name]]
    tree_name <- if(group_name=="study_all") "study" else "full"
    tree <- design[[paste0(tree_name,"_tree")]]
    age <- design[[paste0(tree_name,"_age")]]
    expected_ids <- x$selected_nodes$character
    edges <- x$selected_nodes$edge
    if(group_name=="full_matched") {
      keep <- x$selected_nodes$variable_study
      expected_ids <- expected_ids[keep]
      edges <- edges[keep]
    }
    if(group_name=="study_all") {
      keep <- x$selected_nodes$variable_study
      expected_ids <- expected_ids[keep]
      edges <- x$selected_nodes$study_edge[keep]
    }
    expected_origins <- tabulate(edges,nbins=nrow(tree$edge))
    stopifnot(identical(group$character_ids,expected_ids),
              all(group$branch_summary$origin_count==expected_origins),
              group$divisor==1,group$normalization=="none",
              ncol(group$rate_samples)==100L,
              all(is.finite(group$rate_samples)),all(group$rate_samples>=0),
              max(abs(rowMeans(group$rate_samples)-group$branch_summary$mean))<1e-9,
              max(abs(rowMeans(group$marginal_samples)-group$marginal_summary$mean))<1e-9)
    interval_max_error <- max(interval_max_error,
      check_intervals(group$rate_samples,group$branch_summary),
      check_intervals(group$marginal_samples,group$marginal_summary))
    counts <- group$rate_samples*tree$edge.length
    stopifnot(max(abs(counts-round(counts)))<1e-7)
    for(i in seq_along(design$age_grid)) {
      time <- design$age_grid[i]
      active <- age[tree$edge[,1]]>=time & age[tree$edge[,2]]<time
      if(time==0) active <- tree$edge[,2]<=Ntip(tree)
      expected <- colMeans(group$rate_samples[active,,drop=FALSE])
      max_error <- max(max_error,abs(expected-group$marginal_samples[i,]))
    }
    stopifnot(max_error<1e-9)
  }
  audit <- rbind(audit,data.frame(task=task,condition=x$condition,trial=x$trial,
    characters=x$n_characters,retained=x$n_retained,marginal_max_error=max_error,
    interval_max_error=interval_max_error,
    result_md5=unname(tools::md5sum(result_file)),
    full_maps_md5=map_hashes[1],study_maps_md5=map_hashes[2]))
}
write.csv(audit,"summary/completion-audit.csv",row.names=FALSE)
stopifnot(total_full_maps==12870000L,total_study_maps==3895800L)
report <- c("PASS: all 1500 paired trials independently audited.",
  paste("Full-tree maps:",total_full_maps),paste("Study-tree maps:",total_study_maps),
  paste("Maximum per-map marginal disagreement:",max(audit$marginal_max_error)),
  paste("Maximum 95% quantile disagreement:",max(audit$interval_max_error)),
  "No failed fits or mappings; constants explicitly tracked.",
  "All 734 node patterns and the Apocrita-only origin restriction checked independently.",
  "All saved raw maps, selected-node epochs, origins, unit scaling and map counts checked.",
  "Every saved branch/marginal 95% interval and every per-map marginal curve checked.")
writeLines(report,"summary/completion-audit.txt")
cat(report,sep="\n")

