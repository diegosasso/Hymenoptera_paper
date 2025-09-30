
#--------------------#

### Get ED model for a pair of traits (controlling + dependent character) ###
get_ED_model <- function(mat) {

  # Get character state tokens.
  charst_tk <- apply(mat, 2, function(x) sort(unique(x[x != "?" & x != "-"])) , simplify = F)

  # Get individual models.
  Q_list <- mapply(x = charst_tk, y = 1:2, function(x,y) initQ(x, rate.param = y) , SIMPLIFY = F)

  # Get amalgamated model.
  Q_model <- amaED(Q_list[[1]], Q_list[[2]], type = "ql")

  # Rename character state tokens.
  colnames(Q_model) <- rownames(Q_model) <- 0:(length(unlist(charst_tk)) - 2)

  # Return results.
  return(Q_model)

}

#--------------------#

### Recode non-dependent characters (just to standardize tokens) ###
recode_simple <- function(mat) {
  
  # Get character state tokens.
  charst_tk <- apply(mat, 2, function(x) names(table(x)) )

  # Remove inapplicable and missing tokens.
  old_states <- lapply(charst_tk, function(x) x[!x == "?" & !x == "-"] )

  # Get new states.
  new_states <- lapply(old_states, function(x) paste(0:(length(x) - 1)) )

  # Set token for missing.
  # RevBayes.
  # miss_state <- lapply(new_states, function(x) paste0(x, collapse = " ") )
  # corHMM.
  miss_state <- "?"

  # Recode all characters.
  for (i in 1:dim(mat)[2]) { mat[[i]] <- new_states[[i]][match(mat[[i]], old_states[[i]])] }

  # Recode missings.
  mat[is.na(mat)] <- miss_state
  
  # Return #
  #return(list(matrix = mat, old_states = old.states, new_states = new.states))
  return(mat)

}

#--------------------#

### Recode simple dependencies (binary controlling + dependent character with N states) ###
recode_dep_ql <- function(mat) {

  # Get miss state.
  miss_st <- paste(1:length(unique(mat[[2]][mat[[2]] != "-" & mat[[2]] != "?"])), collapse = "&")

  # Combine characters.
  comb_char <- apply(mat, 1, paste0, collapse = "")

  # Recode state 0.
  comb_char <- gsub(comb_char, pattern = ".\\-", replacement = "0")

  # Recode missings in the first character.
  comb_char <- gsub(comb_char, pattern = "\\?.", replacement = "?")

  # Recode missings in the second character.
  comb_char <- gsub(comb_char, pattern = ".\\?", replacement = miss_st)

  # Set new tokens for character states.
  new_tokens <- unique(comb_char[grep(comb_char, pattern = "\\d{2}")])
  new_tokens <- cbind(new_tokens, 1:length(new_tokens))

  # Recode amalgamated character.
  for (i in 1:dim(new_tokens)[1]) { comb_char[comb_char == new_tokens[i,1]] <- new_tokens[i,2] }

  # Return recoded character.
  new_name <- paste0(names(mat), collapse = "_")

  return(set_names(tibble(X = comb_char), new_name))

}

#--------------------#

### Recode complex dependencies (binary controlling + two dependent binary characters) ###
recode_dep_ql_triple <- function(mat, Q) {

  # Get standard state names.
  st_names <- colnames(Q)[-1]
  st_names <- setNames(st_names, 1:length(st_names))

  # Get miss state.
  miss_st <- paste0(names(st_names), collapse = "&")

  # Combine characters.
  comb_char <- apply(mat, 1, paste0, collapse = "")

  j = 1
  # Recode standard states.
  for (j in 1:length(st_names)) { comb_char <- gsub(comb_char, pattern = st_names[[j]], replacement = j) }

  # Recode state 0.
  comb_char <- gsub(comb_char, pattern = "0..", replacement = "0")
  comb_char <- gsub(comb_char, pattern = ".\\-\\-", replacement = "0")

  # Recode missings in the first character.
  comb_char <- gsub(comb_char, pattern = "\\?..", replacement = "?")

  # Recode missings in the second and third characters.
  comb_char <- gsub(comb_char, pattern = ".\\?\\?", replacement = miss_st)
  comb_char <- gsub(comb_char, pattern = ".\\-\\?", replacement = miss_st)
  comb_char <- gsub(comb_char, pattern = ".\\?\\-", replacement = miss_st)

  # Recode partial inapplicables or missings.
  comb_char <- gsub(comb_char, pattern = "1\\-.", replacement = miss_st)
  comb_char <- gsub(comb_char, pattern = "1\\?.", replacement = miss_st)
  #comb_char <- gsub(comb_char, pattern = "1\\-0", replacement = "1&3")
  #comb_char <- gsub(comb_char, pattern = "1\\?0", replacement = "1&3")
  #comb_char <- gsub(comb_char, pattern = "1\\-1", replacement = "2&4")
  #comb_char <- gsub(comb_char, pattern = "1\\?1", replacement = "2&4")

  comb_char <- gsub(comb_char, pattern = "1.\\-", replacement = miss_st)
  comb_char <- gsub(comb_char, pattern = "1.\\?", replacement = miss_st)
  #comb_char <- gsub(comb_char, pattern = "10\\-", replacement = "1&2")
  #comb_char <- gsub(comb_char, pattern = "10\\?", replacement = "1&2")
  #comb_char <- gsub(comb_char, pattern = "11\\-", replacement = "3&4")
  #comb_char <- gsub(comb_char, pattern = "11\\?", replacement = "3&4")

  # Return recoded character.
  return(set_names(tibble(X = comb_char), paste0(names(mat), collapse = "_")))

}

#--------------------#

## Get number of discrete states ##
get_state_n <- function(char) {
  
  et <- grep(char, pattern = "&")
  char[et] <- "*"
  n_st <- length(unique(char[!char %in% c("?", "-", "*")]))
  
  return(n_st)
  
}

#--------------------#

## Get priors for stochastic mapping ##
get_priors <- function(tree, char, model_obj, model, root_p, dep_char = F) {
  
  # Make state labels.
  st_lb <- paste0(0:(get_state_n(char[,2])-1))
  
  # Get index matrix.
  Q_i <- model_obj$index.mat
  
  # Get Q matrix.
  Q <- model_obj$solution
  
  # Get reconstruction.
  p <- sapply(1:max(Q_i, na.rm = TRUE), function(x) na.omit(c(Q))[na.omit(c(Q_i) == x)][1])
  
  if (dep_char) {

    anc_r <- suppressWarnings(ancRECON_ROOT(phy = tree, data = char, p = p, rate.mat = model, method = "scaled", rate.cat = 1, root.p = root_p, get.likelihood = T, get.tip.states = T))

  } else {
    
    anc_r <- ancRECON_ROOT(phy = tree, data = char, p = p, model = model, method = "scaled", rate.cat = 1, root.p = root_p, get.likelihood = T, get.tip.states = T)
    
  }
  
  # Get root prior.
  if (root_p == "flat") { anc_r$root <- rep((1/dim(Q)),dim(Q)) }
  if (length(anc_r$root) == 0) { root_prior <- anc_r$root } else { root_prior <- setNames(anc_r$root, st_lb) }
  
  # Get tips prior.
  tips_prior <- apply(anc_r$lik.tip.states, 1, function(x) x/sum(x))
  tips_prior <- t(tips_prior)
  rownames(tips_prior) <- tree$tip.label
  colnames(tips_prior) <- st_lb
  
  # Reorganize Q matrix for phytools.
  diag(Q) <- -rowSums(Q, na.rm = T)
  Q[is.na(Q)] <- 0
  colnames(Q) <- rownames(Q) <- st_lb
  
  # Return results.
  return(list("Q_matrix" = Q, "root_prior" = root_prior, "tips_prior" = tips_prior))
  
}

#--------------------#

### Wrapper function for model fitting ###
fit_models <- function(tree, ch, root_p = c("maddfitz", "yang", "flat"), model = NULL, dep_char = F) {
  
  # Extract taxon names.
  taxa <- tree$tip.label
  
  # Get character vector.
  char <- cbind(taxa,ch)
  
  # Set root prior.
  if (root_p == "flat") { root_p_n <- NULL } else { root_p_n <- root_p }
  
  if(dep_char) {
    
    # Fit model.
    best_model <- corHMM(phy = tree, data = char, rate.mat = model, rate.cat = 1, node.states = "none", root.p = root_p_n)
    
    # Get all priors.
    priors <- get_priors(tree = tree, char = char, model_obj = best_model, model = model, root_p = root_p, dep_char = dep_char)
    
  } else {
    
    # Set candidate models.
    models <- c("ER", "SYM", "ARD")
    
    # Set a list to store fitted models.
    fit_mod <- vector(mode = "list", length = length(models))
    
    # Fit models.
    for (j in 1:length(models)) {
      
      fit_mod[[j]] <- corHMM(phy = tree, data = char, model = models[[j]], rate.cat = 1, node.states = "none", root.p = root_p_n)
      #fit_mod[[j]] <- rayDISC(phy = tree, data = char, model = models[[j]], node.states = "none", root.p = root_p_n)
      
    }
    
    # Get best model.
    fits <- aic.w(sapply(fit_mod, function(x) x$AICc))
    best <- min(which(fits == max(fits)))
    best_model <- fit_mod[[best]]
    model <- models[[best]]
    
    print("OK. Passed model fitting.")
    
    # Get all priors.
    priors <- get_priors(tree = tree, char = char, model_obj = best_model, model = model, root_p = root_p, dep_char = dep_char)
    
    print("OK. Passed prior setting.")
  
  }

  # Return results.
  return(c(list("best_model" = best_model), priors))
  
}

#--------------------#

### Wrapper function for quality control of model fitting ###
model_fit_qc <- function(tree, mat, fitted_models, root_p = "maddfitz", thres, trials, dep_char = F, Q_index = NULL) {
  
  # Get putative outliers.
  q_r_high <- which(sapply(fitted_models, function(x) any(rowSums(x$best_model$solution, na.rm = T) > thres) ))
  q_r_high <- names(q_r_high)
  
  if (length(q_r_high) == 0) { return(print("No character detected with rates above the threshold.")) } else {
    
    # Retry model fitting a couple of 'trial' times or until getting better estimates (< 'thres').
    for (i in q_r_high) {
      
      if (dep_char) {
        
        # Model fitting (with corHMM).
        fitted_models[[i]] <- fit_models(tree = hym_tree, ch = mat[[i]], root_p = root_p, dep_char = dep_char, model = Q_index[[i]])
        
        if ( all(rowSums(fitted_models[[i]]$best_model$solution, na.rm = T) < thres) ) { break }
        
      } else {
        
        for (j in 1:trials) {
          
          # Model fitting (with corHMM).
          fitted_models[[i]] <- fit_models(tree = hym_tree, ch = mat[[i]], root_p = root_p)
          
          if ( all(rowSums(fitted_models[[i]]$best_model$solution, na.rm = T) < thres) ) { break }
          
        }
        
      }
      
    }
    
    # Check for outliers again.
    q_r_high_f <- which(sapply(fitted_models, function(x) any(rowSums(x$best_model$solution, na.rm = T) > thres) ))
    q_r_high_f <- names(q_r_high_f)
    
    # Return results.
    if(length(q_r_high_f) < length(q_r_high)) { 
      
      print("All good! Some issues solved.")
      return(list("fitted_models" = fitted_models, "outliers" = q_r_high_f)) 
      
    } else {
      
      print("Issues persist. These characters might be outliers.")  
      return(list("fitted_models" = fitted_models, "outliers" = q_r_high_f))
      
    }
    
  }
  
}

#--------------------#

### Wrapper function for processing and recoding dependent characters ###
process_dep_data <- function(mat, char_annot, complex_dep = F) {
  
  # Extract controlling characters.
  crtl_chars <- char_annot %>% filter(dep_tag == "contrchar")
  
  # Extract dependency groups and organize data.
  dep_groups <- char_annot %>% filter(dep_tag != "contrchar") %>% group_by(depends_on) %>% group_split()
  names(dep_groups) <- sapply(dep_groups, function(x) unique(x$depends_on) )
  dep_groups <- dep_groups[order(as.numeric(gsub(names(dep_groups), pattern = "CH", replacement = "")))]
  
  # Get dependency group names.
  char_names <- names(dep_groups)
  
  # Create a list to store information on dependent characters.
  dep_data_info <- vector(mode = "list", length = length(dep_groups))
  
  # Create a list to store merged characters.
  mat_merg <- vector(mode = "list", length = length(dep_groups))
  
  # Create a list to store Q matrices indices.
  Q_list <- vector(mode = "list", length = length(dep_groups))
  
  i = 1
  # Loop over all dependencies.
  for (i in 1:length(dep_groups)) {
    
    cat(paste0("\n", "Working on group: ", char_names[[i]], ": ", Sys.time(), "\n"))
    
    # Extract data on dependency.
    char_dat <- dep_groups[[char_names[[i]]]]
    
    # Get dependent group.
    dep_group <- c(unique(char_dat$depends_on), char_dat$char_id)
    
    # Process dependent character info.
    x <- char_annot %>% filter(char_id %in% dep_group) %>% select(char_id, onto_id, label, type)
    x <- apply(x, 2, function(y) paste0(unique(y), collapse = "_") )
    
    # Store dependent character info.
    dep_data_info[[i]] <- x
    
    # Get character vectors.
    char_vecs <- mat[dep_group]
    
    # Check if trigger state is 0 or 1 in the controlling character. Invert coding if necessary.
    if (unique(char_dat$upon_state) == 0) {
      
      char_vecs[[1]] <- gsub(char_vecs[[1]], pattern = "0", replacement = "X")
      char_vecs[[1]] <- gsub(char_vecs[[1]], pattern = "1", replacement = "0")
      char_vecs[[1]] <- gsub(char_vecs[[1]], pattern = "X", replacement = "1")
      
    }
    
    if (complex_dep) {
      
      # PROCESS COMPLEX DEPENDENCIES #
      
      # Get possible state combinations in the dependent characters.
      char_comb <- apply(char_vecs, 1, paste0, collapse = "")
      char_comb <- unique(char_comb[grep(char_comb, pattern = "\\d{3}")])
      
      # Set standard model.
      Q <- initQ(c("0", "1"), c(1,2))
      Q <- amaED(Q, amaSMM(Q,Q), type = "ql")
      
      # Check if all states are observed. If not then get a reduced model.
      if (length(char_comb) != 4) {
        
        # Set reduced model.
        Q <- Q[colnames(Q) %in% c("0", char_comb),colnames(Q) %in% c("0", char_comb)]
        
      }
      
      # Merge and recode dependent characters. Store recoded characters.
      mat_merg[[i]] <- recode_dep_ql_triple(char_vecs, Q)
      
      # Store index matrices.
      colnames(Q) <- rownames(Q) <- paste0(0:(length(colnames(Q))-1))
      diag(Q) <- 0
      Q[is.na(Q)] <- 0
      Q_list[[i]] <- Q
      
    } else {
      
      # PROCESS SIMPLE DEPENDENCIES #
      
      # Merge and recode dependent characters. Store recoded characters.
      mat_merg[[i]] <- recode_dep_ql(char_vecs)
      
      # Get amalgamated Q matrix indices. Store index matrices.
      Q <- get_ED_model(char_vecs)
      diag(Q) <- 0
      Q[is.na(Q)] <- 0
      Q_list[[i]] <- Q
      
    }
    
  }
  
  # Build tibble with dependent character info.
  dep_data_info <- do.call(bind_rows, dep_data_info)
  
  # Build new matrix with merged characters.
  mat_merg <- do.call(bind_cols, mat_merg)
  
  # Rename list with fitted models.
  names(Q_list) <- dep_data_info$char_id
  
  # Return results.
  return(list("matrix_merged" = mat_merg, "Q_index" = Q_list, "dep_char_info" = dep_data_info))
  
}

#--------------------#

### Wrapper function for stochastic character mapping ###
paramo_stm <- function(tree, fitted_models, start = 1, n_stm = 100, outdir = "stmaps") {
  
  for (i in start:length(fitted_models)) {
    
    tryCatch(
      
      {
        
        R.utils::withTimeout(
          
          {
            
            cat(paste0("\n", "Working on: ", names(fitted_models)[[i]], ": ", Sys.time(), "\n"))
            
            # Get Q matrix and priors.
            Q <- fitted_models[[i]]$Q_matrix
            pi <- fitted_models[[i]]$root_prior
            char <- fitted_models[[i]]$tips_prior
            
            # Sample stochastic maps.
            stm <- suppressWarnings(make.simmap(tree = tree, x = char, Q = Q, pi = pi, nsim = n_stm))
            
            # Save RDS files.
            saveRDS(stm, file = paste0(outdir, "/", names(fitted_models)[[i]], ".RDS"))
            
          }, timeout = 300, onTimeout = "error"
          
        )
        
      }, TimeoutException = function(ex) { message(paste0("Skipped ", names(fitted_models)[[i]], " due to timeout"))
        NULL
      }
      
    )
    
  }
  
}

#--------------------#

### Wrapper function for rate estimation with ontophylo ###
pnhpp_rates <- function(stm, res, b_width) {
  
  # Discretize a reference tree.
  tree_discr <- discr_Simmap(stm[[1]], res = res)
  
  # Merge state categories across branches.
  cat(paste0("\n", "Starting merging state categories: ", Sys.time(), "\n"))
  stm_merg <- merge_tree_cat_list(stm_amalg)
  
  # Calculate Hamming distances.
  cat(paste0("\n", "Starting calculating hamming distances: ", Sys.time(), "\n"))
  path_hm <- path_hamming_over_trees_KDE(stm_merg)
  
  # Make path data.
  path_data <- make_data_NHPP_KDE_Markov_kernel(path_hm, add.psd = F)
  
  # Copy original path data.
  path_data_c <- path_data
  
  # Estimate bandwidth.
  bdw <- estimate_band_W_mod(tree_discr, path_data_c, band.width = b_width)
  bdw <- mean(bdw)
  
  # Kernel Density Estimator (KDE).
  cat(paste0("\n", "Starting estimating KDEs: ", Sys.time(), "\n"))
  edge_KDE <- estimate_edge_KDE(tree_discr, Path.data = path_data_c, h = bdw)
  
  # Smoothing and normalization of KDE data.
  edge_KDE$Maps.mean.loess <- suppressWarnings(loess_smoothing_KDE(tree_discr, edge_KDE))
  edge_KDE$Maps.mean.loess.norm <- normalize_KDE(tree_discr, edge_KDE$Maps.mean.loess)
  
  # Calculate the lambda statistics.
  lambda_post <- posterior_lambda_KDE(stm_merg)
  
  # Approximate posterior distribution.
  edge_KDE$lambda.mean <- make_postPois_KDE(edge_KDE$Maps.mean.norm, lambda_post, lambda.post.stat = "Mean")
  edge_KDE$lambda.mean.loess <- make_postPois_KDE(edge_KDE$Maps.mean.loess.norm, lambda_post, lambda.post.stat = "Mean")
  
  # Make data for contmaps.
  nhpp_lambda_mean <- make_contMap_KDE(tree_discr, edge_KDE$lambda.mean.loess)
  
  # Make data for edge profiles.
  edge_profs_lambda_mean <- edge_profiles4plotting(tree_discr, edge_KDE$lambda.mean.loess)
  
  # Return results.
  return(list("stm" = stm_merg, "hamming" = path_hm, "path_data" = path_data,
              "KDE" = edge_KDE, "contmap" = nhpp_lambda_mean,
              "edgeplot" = edge_profs_lambda_mean))
  
}

#--------------------#

# Wrapper function for making morphospace data.
make_morphospace <- function(stm, n_samples) {
  
  # Merge state categories across branches.
  stm_merg <- merge_tree_cat_list(stm)
  
  # Get a sample of stochastic maps.
  stm_merg <- stm_merg[seq(from = (length(stm_merg)/n_samples), to = length(stm_merg), by = n_samples)]
  
  # Multidimensional scaling.
  MDS_list <- lapply(stm_merg, MultiScale.simmap)
  
  # Return results.
  return(MDS_list)
  
}

#--------------------#

### Wrapper function for making edgeplots ###
make_edgeplot <- function(cont_data, edg_data, tip_lb, ana_name) {
  
  # Workaround for improve plotting of zero-rates.
  if (any(edge_data$Y <= 0)) { edge_data$Y <- replace_zero(edge_data$Y) }
  
  # Get tree height.
  Tmax <- max(nodeHeights(cont_data$tree))
  
  # Set plot layout.
  layout(matrix(c(1,2),ncol = 1), heights = c(2,1))
  
  # Plot contmap.
  plot.contMap(cont_data, lwd = 3, outline = F, legend = F, ftype = "off", plot = F, mar = c(0.1, 3.45, 0.1, 0.35))
  
  # Add tip labels.
  tiplabels(pch = 19, col = tip_lb, cex = 0.6)
  
  # Plot edgeplot.
  plot_edgeprof <-
    
    ggplot(data = edge_data, aes(x = X-Tmax, y =  Y, group = edge.id, color = Y)) +
    
    geom_line(alpha = 1, linewidth = 0.5) +
    
    scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = 0.7)) ) +
    
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16),
          plot.margin = unit(c(2.3,0.87,0.1,0.1), 'cm'),
          legend.position = 'none') +
    
    xlab('time') + ylab('rate') +
    
    scale_x_continuous(limits = c(-round(Tmax + 5, 0), 0),
                       breaks = -1*seq(from = 0, to = Tmax, by = Tmax/5) %>% round(0),
                       labels = seq(from = 0, to = Tmax, by = Tmax/5) %>% round(0) ) +
    scale_y_continuous(limits = c(0, round(max(nhpp$edgeplot$Y)*1.2, 3))) +
    coord_cartesian(expand = FALSE)
  
  vp <- grid::viewport(height = unit(0.5,"npc"), width = unit(1, "npc"), just = c("left",'top'), y = 0.5, x = 0)
  
  print(plot_edgeprof, vp = vp)
  
  title(main = paste0(ana_name), font.main = 2, line = -0.5, cex.main = 0.8)
  
}

#--------------------#

### Wrapper function for making morphospace plots ###
make_morphospace_plot <- function(MDS, ana_type = c("groups", "time"), tip_tags, 
                                  add_noise = F, noise = c(0.3, 0.3)) {
  
  if (add_noise == TRUE) { MDS <- add_noise_MD(MDS, add.noise = noise) }
  
  # Get tip ids.
  tip_ids <- which(MDS$Points$sp_extant == "yes")
  
  # Update tibbles.
  # Add tip ids.
  MDS$Points <- mutate(MDS$Points, tip.id = c(1:nrow(MDS$Points)))
  
  # Add tip and non-tip tags
  MDS$Points <- MDS$Points %>% mutate(groups = rep("no_tip", dim(MDS$Points)[1]))
  MDS$Points$groups[tip_ids] <- tip_tags
  
  # Set some parameters for the plots.
  # Tree height.
  Tmax = max(MDS$Points$time)
  
  #####################
  ### MAKE MDS PLOT ###
  #####################
  
  # Base plot. Edges.
  base_plot <- geom_segment(data = MDS$Lines, aes(x = start.V1, y = start.V2, xend = end.V1, yend = end.V2), 
                            colour = "red", linewidth = 0.3, linetype = 1, alpha = 0.3)
  
  # Add points. Ancestral states.
  points_anc_plot <- geom_point(data = MDS$Points %>% filter(sp_extant == "no"), 
                                aes(x = V1, y = V2), color = "grey", 
                                alpha = 0.5, size = 1.5, show.legend = TRUE)
  
  # Add points. Present states. Color by group.
  points_pres_plot <- geom_point(data = MDS$Points %>% filter(sp_extant == "yes"), 
                                 aes(x = V1, y = V2, color = groups), alpha = 0.8, size = 3)
  
  # Add points. Ancestral and present states. Color by time.
  points_time_plot <- geom_point(data = MDS$Points, aes(x = V1, y = V2, color = time), alpha = 0.8, size = 3)
  
  if (ana_type == "groups") {
    
    # MDS plot. By group.
    MDS_plot <- ggplot() + base_plot + points_anc_plot + points_pres_plot + 
      
      xlab("Coor1") + ylab("Coor2") +
      
      cowplot::theme_cowplot() + 
      
      theme(legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    
  }
  
  if (ana_type == "time") {
    
    # MDS plot. By time.
    MDS_plot <- ggplot() + base_plot + points_time_plot + 
      
      scale_color_gradient(low = "white", high = "blue", limits = c(0, Tmax), 
                           breaks = seq(0, Tmax, 50), labels = seq(0, Tmax, 50) %>% rev) + 
      
      xlab("Coor1") + ylab("Coor2") +
      
      cowplot::theme_cowplot() + 
      
      theme(legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    
  }
  
  # Return plot.
  return(MDS_plot)
  
}

#--------------------#

# internal function: replace zero or near-zero rate values.
replace_zero <- function(x) { 
  
  x[x <= 0] <- 1e-10
  
  return (x) 
  
}

#--------------------#

##################################################
# PATCHED FUNCTIONS FROM ONTOPHYLO
##################################################


estimate_band_W_mod <- function (tree.discr, data.path, band.width = c("bw.nrd0", "bw.nrd0", "bw.ucv", "bw.bcv", "bw.SJ")) {
  
  otus <- c(1:length(tree.discr$tip.label))
  Edges <- match(otus, tree.discr$edge[, 2])
  Band.width <- c()
  
  for (i in Edges) {
    
    # Skip a path if not enough changes were observed (update).
    if (length(data.path[[i]]) <= 3 ) { next } else { dt <- data.path[[i]] }
    
    if (band.width == "bw.nrd") 
      h <- bw.nrd(dt)
    if (band.width == "bw.nrd0") 
      h <- bw.nrd0(dt)
    if (band.width == "bw.ucv") 
      h <- bw.ucv(dt, lower = 0.01, upper = 20)
    if (band.width == "bw.bcv") 
      h <- bw.bcv(dt, lower = 0.01, upper = 20)
    if (band.width == "bw.SJ") 
      h <- bw.SJ(dt, lower = 0.01, upper = 20)
    
    Band.width <- c(Band.width, h)
    
  }
  
  return(Band.width)

}
