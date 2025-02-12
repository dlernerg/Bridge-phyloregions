
sim_out_fun <- function(i,bridge){
  #sim_num <- sim_num
  stop_max_spec <- stop_max_spec
  min_spec <- 20
  N.eqR <- 1
  max.sim <- max.sim
  sim_out_list_par <- list()
  
  sim_fun_bridge <- function(sim_vec,stop_max_spec,min_spec){
    require(dplyr)
    require(phytools)
    require(ggplot2)
    require(picante)
    
    # functions ---------------------------------------------------------------
    initialize_species_matrix <- function(num_regions) {
      matrix(ncol = num_regions, nrow = 1, data = c(rep(1,num_regions)))
    }
    
    disp_xtc.fun <- function(species_vector) {
      # Initialize disp.ext_event matrix for multiple regions
      disp.ext_event <- matrix(0, nrow = 3, ncol = num_regions)
      
      # Loop through each region to calculate dispersal and extinction events
      for (region in 1:num_regions) {
        #if moved to bridge region then identify where from
        
        if (species_vector[region] == 1) {
          # Determine adjacent regions
          prev_region <- region - 1
          next_region <- region + 1
          
          # Dispersal to previous region (if it's empty and exists)
          if (prev_region >= 1 && species_vector[prev_region] == 0) {
            disp.ext_event[1, prev_region] <- rbinom(1, 1, disp_prob_matrix[region, prev_region])
            #this will allow to identify movement to the bridge
            if (prev_region == 2 && disp.ext_event[1,prev_region] == 1){
              disp.ext_event[3,region] <- 1
            }
          }
          
          # Dispersal to next region (if it's empty and exists)
          if (next_region <= num_regions && species_vector[next_region] == 0) {
            disp.ext_event[1, next_region] <- rbinom(1, 1, disp_prob_matrix[region, next_region])
            if (next_region == 2 && disp.ext_event[1,next_region] == 1){
              disp.ext_event[3,region] <- 1
            }
          }
          
          disp.ext_event[2, region] <- -rbinom(1, 1, prob_xtc[region])
          
        }
      }
      return(disp.ext_event)
    }
    
    
    
    # update_species_matrix <- function(all_species_matrix, disp_ext_event) {
    #   # Assuming all_species_matrix and disp_ext_event are matrices with the same number of columns (regions)
    #   num_regions <- ncol(all_species_matrix)
    #   
    #   # Updating the species matrix with dispersal and extinction events
    #   for (region in 1:num_regions) {
    #     # Update the species count for each region
    #     # Adding dispersal events and subtracting extinction events
    #     all_species_matrix[region] <- max(0, all_species_matrix[region] + sum(disp_ext_event[1, ]) - disp_ext_event[2, region])
    #   }
    #   
    #   return(all_species_matrix)
    # }
    # 
    find_speciation_events <- function(species_vector) {
      speciation_events <- rep(FALSE, num_regions)
      
      # Checking each region for the presence of species and potential speciation
      for (region in 1:num_regions) {
        if (species_vector[region] > 0) {
          # Assume some logic or condition to determine if speciation occurs
          # For example, based on a probability or a certain condition
          speciation_events[region] <- rbinom(1,1,prob_spec[region]) # Example condition
        }
      }
      
      return(speciation_events)
    }
    
    prepare.output_fun <- function(spec.id,next_gen){
      
      spec.id$spec_event <- next_gen$n
      
      from.node <- spec.id$species_origin
      to.node <- spec.id$ID
      edge <- data.frame(cbind(from.node,to.node),spec_event = spec.id$spec_event,spp = spec.id$ID, alive = spec.id$alive)
      
      edge.length <- c()
      for (i in 1:(nrow(spec.id))){
        spec.origin <- spec.id$species_origin[i]
        edge.length[i] <- spec.id$generation[i]-spec.id$generation[which(spec.id$ID == spec.origin)]
      }
      
      stem.depth <- spec.id$generation
      
      edge$spp <- NA
      node <- 1
      spp <- 1
      
      a <- which(!is.na(all_species_matrix_f[,1]))
      b <- which(!is.na(all_species_matrix_f[,2]))
      d <- which(!is.na(all_species_matrix_f[,3]))
      a <- unique(c(a,b,d))
      edge$alive <- 0
      edge$alive[a] <- 1
      
      for (i in 1:nrow(edge)){
        
        if(edge$alive[i]==1){
          edge$node[i] = NA
          edge$spp[i] <- spp 
          spp <- spp + 1
          
        }else{
          edge$node[i] = node 
          node <- node + 1
        }
      }
      
      edge.f <- edge %>% 
        filter(spec_event > -999 | alive == 1) %>%
        mutate(ID = to.node)
      
      run<- unique(edge.f$from.node)
      
      n = -1
      for (i in run) {
        if (any(edge.f[, 'from.node'] == i)) {
          edge.f[which(edge.f[, 'from.node'] == i), 'from.node'] = n
          edge.f[which(edge.f[, 'to.node'] == i), 'to.node'] = n
          n = n - 1
        }
      }
      
      edge.f.f <- edge.f %>% 
        filter(!(to.node > 0 & alive == 0))
      
      
      filtered <- which(!(edge$to.node %in% edge.f.f$ID))
      
      if(identical(filtered,integer(0))){
        edge.out <- edge
        edge.length.out <- edge.length
        edge.length.out[1] <- 0
        stem.depth.out <- stem.depth
      }else{
        edge.out <- edge[-filtered,]
        edge.length.out <- edge.length[-filtered]
        edge.length.out[1] <- 0
        stem.depth.out <- stem.depth[-filtered]  
      }
      
      return(list(edge = edge.out,edge.length = edge.length.out,stem.depth = stem.depth.out))
    }
    
    make_phylo_fun <- function(output.list){
      edge.out <- output.list$edge
      edge.length.out <- output.list$edge.length
      stem.depth.out <- output.list$stem.depth
      
      edge.length.out[which(edge.out$alive == 1)] = sim_num.new - stem.depth.out[which(edge.out$alive == 1)]
      run<- unique(edge.out$from.node)
      edge.out.m = as.matrix(edge.out)
      n = -1
      
      for (i in run) {
        if (any(edge.out.m[, 'from.node'] == i)) {
          edge.out.m[which(edge.out.m[, 'from.node'] == i), 'from.node'] = n
          edge.out.m[which(edge.out.m[, 'to.node'] == i), 'to.node'] = n
          n = n - 1
        }
      }
      
      edge.out.m[which(edge.out.m[, c('from.node', 'to.node')] > 0)] = 1:length(edge.out.m[which(edge.out.m[, c('from.node', 'to.node')] > 0)])
      tip.label = edge.out.m[, 'to.node'][!is.na(edge.out[, 'spp'])]
      edge.only = edge.out.m[, c('from.node', 'to.node')]
      mode(edge.only) = "character"
      mode(tip.label) = "character"
      obj = list(edge = edge.only, edge.length = edge.length.out, tip.label = tip.label)
      class(obj) = "phylo"
      phylo.out = old2new.phylo(obj)
      return(phylo.out)
    }
    
    change_matrix_dimensions <- function(mat, nrow) {
      mat <- mat[1:nrow, 1:dim(mat)[2]]
      return(mat)
    }
    
    exp_seq.fun <- function(start,end,length){
      return(exp(seq(log(start), log(end), length.out = length)))
    }
    
    
    
    
    # load variables --------------------------------------------------------------
    
    
    stop_max_spec <- stop_max_spec
    min_spec <- min_spec
    
    num_regions <- 3
    
    prob_spec <- c(sim_vec$spec_trop,NA,sim_vec$spec_temp)
    prob_xtc <- c(sim_vec$xtc_trop,NA,sim_vec$xtc_temp)
    disp_prob_matrix <- matrix(ncol = num_regions, nrow = num_regions,NA)
    disp_prob_matrix[1,2] <- sim_vec$disp_trop
    disp_prob_matrix[2,1] <- sim_vec$disp_btrop
    disp_prob_matrix[2,3] <- sim_vec$disp_btemp
    disp_prob_matrix[3,2] <- sim_vec$disp_temp
    
    stop_max_spec <- stop_max_spec
    min_spec <- min_spec
    
    
    # run sim -----------------------------------------------------------------
    
    species_vector <- initialize_species_matrix(num_regions)
    species_vector[1,2] <- 0
    
    alive_species_count_trop <- c()
    alive_species_count_bridge <- c()
    alive_species_count_temp <- c()
    
    
    all_species_matrix <- matrix(NA, nrow = empty_row, ncol = num_regions)
    all_species_matrix[1,] <- species_vector
    #all_species_matrix <- list()
    #all_species_matrix[[1]] <- species_vector
    spec.id <- data.frame(species_origin = 0,
                          generation = 0,
                          source = 1,
                          ID = 1, 
                          alive = 1, 
                          speciated = NA
    )
    
    next_gen <- data.frame(n=-999,spec = 1)
    
    sim.num <- 1
    N.eqR <- N.eqR
    disp_fun_triggered <- FALSE
    alive_species_count <- c()
    prob_evo <- list(c(NA,NA))
    #prob_evo[[1]] <- c(prob_spec[1],prob_xtc[1])
    
    repeat{
      
      n_spec <- nrow(spec.id) 
      new.spec <- 1
      alive_sum <- sum(spec.id$alive)
      
      use_disp_fun <- alive_sum >= N.eqR
      
      if (alive_sum > stop_max_spec || alive_sum == 0 || sim.num > max.sim) {
        break
      } 
      
      for (j in 1:n_spec){
        
        #make sure it makes sense to put disp_to_bridge here
        disp_to_bridge = FALSE
        
        species_vector <- all_species_matrix[j,]
        
        prob_evo_at_j <- c()
        
        if (isTRUE(all(is.na(species_vector))) || sum(species_vector) == 0) {
          #all_species_matrix[j,] <- matrix(ncol = 2, nrow = 1, c(NA, NA))
          # spec.id$alive[j] <- 0
          # spec.id$status[j] <- "dead"
          next 
        }
        
        #load the probabilities of extinction and speciation for the bridge region of this species.
        if (is.null(prob_evo[[j]])){
          prob_spec[2] <- NA
          prob_xtc[2] <- NA
        } else {
          prob_spec[2] <- prob_evo[[j]][1]
          prob_xtc[2] <- prob_evo[[j]][2]
          
        }
        
        
        ##this first step define the dispersal extinction and speciation events
        
        #chance of a dispersal or extinction event
        if (use_disp_fun || disp_fun_triggered) {
          disp.ext_event <- disp_xtc.fun(species_vector)
          
          #the third row in the disp.ext_event is required for the identification of where the disp to the bridge region came from 
          if (any(disp.ext_event[3,1:3] == 1)){
            disp_to_bridge = TRUE
            #i assigned the [1] for the disp_from given the unique and uncommon case that the dispersal to bridge will come from both trop and temp 
            disp_from <- which(disp.ext_event[3,] ==1)
            if (length(disp_from>1)){
              disp_from <- disp_from[sample(length(disp_from),1)]
            }
            #if (disp_from == 2){disp_from = 3}
            prob_evo_at_j <- c(prob_spec = prob_spec[disp_from],prob_xtc = prob_xtc[disp_from])
            
            prob_spec[2] <- prob_evo_at_j[1]
            prob_xtc[2] <- prob_evo_at_j[2]
            #assign the prob_spec for this generation for this species, given that we haven't carried out speciation at this region yet
            #prob_spec[2] <- prob_spec[disp_from]
          }
          
          species_vector_new <- species_vector + disp.ext_event[1,] + disp.ext_event[2,]
          
          # Ensure this block always gets executed in future iterations
          disp_fun_triggered <- TRUE
        } 
        else {
          species_vector_new <- species_vector
          
        }
        
        if (sum(species_vector_new) ==0){
          all_species_matrix[j,] <- species_vector_new
          spec.id$alive[j] <- 0
          next
        }
        
        
        speciation_event <- find_speciation_events(species_vector_new)
        
        
        if (any(speciation_event != 0)) {
          species_vector_new_matrix <- matrix(species_vector_new, nrow = 1, ncol = 3)
          
          all_species_matrix[j,] <- matrix(ncol = num_regions, nrow = 1, data = NA)
          spec.id$alive[j] <- 0
          spec.id$speciated[j] = "Y"
          
          speciation_event_time <- speciation_event * round(runif(length(speciation_event)),5)
          
          new_species <- list()  # Temporarily store new species
          
          # This replaces multiple conditions checking speciation_event values
          new_species_counter <- 0
          
          if (sum(speciation_event) == 1){
            
            idx <- which(speciation_event == 1)
            
            new_species_matrix <- matrix(ncol = num_regions, nrow = 1, data = 0)
            new_species_matrix[1, idx] <- 1
            
            new_species <- append(new_species, list(species_vector_new_matrix, new_species_matrix))
            
            new_spec <- data.frame(species_origin = rep(j,2),
                                   generation = rep(sim.num+(speciation_event_time[idx]),2),
                                   source = idx,
                                   ID = n_spec + new.spec + new_species_counter + (0:1),
                                   alive = rep(1,2),
                                   speciated = NA)
            
            spec.id <- bind_rows(spec.id,new_spec)
            
            
            next_gen[j,1] <- next_gen[j,1]+2
            
            new_gen <- data.frame(n = rep(-999,2), spec = n_spec + new.spec + new_species_counter + (0:1))
            
            next_gen <- bind_rows(next_gen,new_gen)
            
            new_species_counter <- new_species_counter + 2
            
          }else {
            
            sorted_times <- order(speciation_event_time)
            spec.id_prev <- spec.id[FALSE,]
            
            
            for (idx in sorted_times) {
              
              if (speciation_event[idx] != 0) {
                spec_origin_in_loop <- j
                
                
                #this is to remove the species that speciates twice
                if (new_species_counter>0){
                  
                  new_species[[new_species_counter-1]] <-  matrix(ncol = num_regions, nrow = 1, data = NA)
                  spec.id_prev$alive[new_species_counter-1] <- 0
                  spec.id_prev$speciated[new_species_counter-1] <- "Y"
                  spec_origin_in_loop <-  n_spec + new.spec + new_species_counter - 2
                  
                }
                
                new_species_matrix <- matrix(ncol = num_regions, nrow = 1, data = 0)
                new_species_matrix[1, idx] <- 1
                
                new_species <- append(new_species, list(species_vector_new_matrix,new_species_matrix))
                #new_species <- append(new_species, list(matrix(ncol = 2, nrow = 1, c(ifelse(idx==1,1,0), ifelse(idx==2,1,0)))))
                
                new_spec <- data.frame(species_origin = rep(spec_origin_in_loop,2),
                                       generation = rep(sim.num+(speciation_event_time[idx]),2),
                                       source = idx,
                                       ID = n_spec + new.spec + new_species_counter + (0:1),
                                       alive = rep(1,2),
                                       speciated = NA)
                
                
                
                spec.id_prev <- bind_rows(spec.id_prev,new_spec)
                
                next_gen[spec_origin_in_loop,1] <- next_gen[spec_origin_in_loop,1]+2
                
                new_gen <- data.frame(n = rep(-999,2), spec = n_spec + new.spec + new_species_counter + (0:1))
                
                next_gen <- bind_rows(next_gen,new_gen)
                
                new_species_counter <- new_species_counter + 2
                
              }
              
            }
            
            spec.id <- bind_rows(spec.id,spec.id_prev)
            
            
          }
          
          
          new.spec <- new.spec+ (sum(speciation_event) * 2)
          new_species <- do.call(rbind, lapply(new_species, function(x) as.matrix(x)))
          
          start <- nrow(spec.id)-(new_species_counter-1)
          end <- nrow(spec.id)
          
          all_species_matrix[start:end,] <- new_species 
          
          #prob_evo[[start:end]] <- prob_evo_at_j
          prob_evo[start:end] <- rep(list(c(NA,NA)), length(end - start + 1))
          
          is_one <- all_species_matrix[start:end, 2] == 1
          is_one[is.na(is_one)] <- F 
          
          if (disp_to_bridge == TRUE){
            
            prob_evo[which(is_one) + start -1] <- rep(list(c(prob_evo_at_j)), length(end - start + 1))
            
          } else if (any(is_one)){
            
            prob_evo[which(is_one) + start -1] <- rep(list(c(prob_evo[[j]])), length(end - start + 1))
            
          }
          
        } else {
          
          all_species_matrix[j,] <- species_vector_new
          #prob_evo[[j]] <- prob_evo_at_j
          
          if (disp_to_bridge == TRUE){
            prob_evo[[j]] <- prob_evo_at_j
          } 
          
        }  
        #now depending if there were or weren't speciation events:
        
        #no speciation events, keep the original matrix, and add a new generation 
        
        
      } 
      
      sim.num<- sim.num+1
      # alive_count_for_current_gen_trop <- sum(spec.id$alive == 1 & spec.id$source == 1)
      # alive_species_count_trop <- c(alive_species_count_trop, alive_count_for_current_gen_trop)
      # 
      # alive_count_for_current_gen_bridge <- sum(spec.id$alive == 1 & spec.id$source == 2)
      # alive_species_count_bridge <- c(alive_species_count_bridge, alive_count_for_current_gen_bridge)
      # 
      # alive_count_for_current_gen_temp <- sum(spec.id$alive == 1 & spec.id$source == 3)
      # alive_species_count_temp <- c(alive_species_count_temp, alive_count_for_current_gen_temp)
      # 
    }
    
    sim_num.new <- sim.num
    col.plot <- palette()[c(4,2)]
    
    if (nrow(spec.id)==1){
      
      #all_species_array =array(unlist(all_species_matrix),dim=c(nrow(all_species_matrix[[1]]),2,length(all_species_matrix)))
      #sum.arr <- apply(all_species_array, 1:2, sum,na.rm = T) %>% as.data.frame() %>% rename(Temperate =  V2, Tropic = V1) %>% reshape2::melt()
      #sum.arr$gen <- rep(1:(nrow(sum.arr)/2),2)
      
      return(list(out = "F", reason = "no species"))
      
      
    }else if (sum(spec.id$alive)<min_spec){
      
      #all_species_array =array(unlist(all_species_matrix),dim=c(nrow(all_species_matrix[[1]]),2,length(all_species_matrix)))
      #sum.arr <- apply(all_species_array, 1:2, sum,na.rm = T) %>% as.data.frame() %>% rename(Temperate =  V2, Tropic = V1) %>% reshape2::melt()
      #sum.arr$gen <- rep(1:(nrow(sum.arr)/2),2)
      
      
      return(list(out = "F",reason = "too few spec", spec.id = spec.id, sim.num = sim.num, total_spec = nrow(spec.id)))
      
    }else if (sum(spec.id$alive)>max_spec){
      
      return(list(out = "F",reason = "too many spec", spec.id = spec.id, sim.num = sim.num, total_spec = nrow(spec.id)))  
      
    }else{
      
      all_species_matrix_f <- all_species_matrix[1:nrow(spec.id),]
      
      #all_species_array =array(unlist(all_species_matrix),dim=c(1,num_regions,length(all_species_matrix)))
      
      #in case you do need to use sum.arr, you will need to find a better way to do so, probably including the count of alive every loop (good way to keep count of the evolution...
      
      #sum.arr <- apply(all_species_array, 1:2, sum,na.rm = T) %>% as.data.frame() %>% setNames(c("Tropic","Temperate")) %>% reshape2::melt()
      #sum.arr$gen <- rep(1:(nrow(sum.arr)/2),2)
      
      
      
      
      # prepare output ----------------------------------------------------------
      
      
      output.list <- prepare.output_fun(spec.id,next_gen)
      
      # make phylo --------------------------------------------------------------
      
      
      
      phylo.sim <- make_phylo_fun(output.list)
      
      #if the phylogenetic tree doesn't work for some reason. 
      continue <- tryCatch({write.tree(phylo.sim)}, error = function(e){    
        
        continue <- "NO"
        return(continue)})
      
      if (continue == "NO"){
        
        return(list(out = "F", reason = "failed phylo"))
        
      }else{
        phylo.sim = read.tree(text = write.tree(phylo.sim))
        
        edge <- output.list$edge
        spp.alive <- data.frame(spp = edge$spp[(edge$alive==1)], ID = edge$to.node[edge$alive==1])
        
        alive_geo <- all_species_matrix_f[spp.alive$ID,]
        rownames(alive_geo) <- spp.alive$spp
        colnames(alive_geo) <- c("Trop","Bridge","Temp") 
        
        #pies <- make.pies(alive_geo, phylo.sim)
        
        phydist <- cophenetic.phylo(phylo.sim)
        mpd <-  ses.mpd(t(alive_geo),phydist, null.model = "taxa.labels" ,runs = 50, abundance.weighted = F)
        mpd.df <- c(MPD.z.1 = mpd$mpd.obs.z[1], MPD.z.2 = mpd$mpd.obs.z[2],MPD.z.3 = mpd$mpd.obs.z[3])
        
        return(list(out = "Y", sim_num = sim_num.new, spec.id = spec.id, next_gen = next_gen, all_species_array = all_species_matrix_f, mpd.df = mpd.df))
        #return(list(out = "Y", sim_num = sim_num.new, tree = phylo.sim, output.list = output.list, spec.id = spec.id, next_gen = next_gen, alive_geo = alive_geo,mpd.df = mpd.df)) 
      }
    } 
  }
  sim_fun_nobridge <- function(sim_vec,stop_max_spec,min_spec){
    require(dplyr)
    require(phytools)
    require(ggplot2)
    require(picante)
    # Functions ---------------------------------------------------------------
    
    disp_xtc.fun <- function(species_vector){
      
      # Initialize disp.ext_event
      disp.ext_event <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("disp", "xtc"), c("Trop", "Temp")))
      
      # Identify current position
      position.current <- which(species_vector == 1)
      
      # Calculate events based on position
      if (length(position.current) == 1) {
        if(position.current == 1){
          disp.ext_event[1, 2] <- rbinom(1, 1, p_disp_trop)
          disp.ext_event[2, 1] <- -rbinom(1, 1, p_xtc_trop)
        } else {
          disp.ext_event[1, 1] <- rbinom(1, 1, p_disp_temp)
          disp.ext_event[2, 2] <- -rbinom(1, 1, p_xtc_temp)
        }
      } else {
        disp.ext_event[2, ] <- -rbinom(2, 1, c(p_xtc_trop, p_xtc_temp))
      }
      
      return(disp.ext_event)    
    }
    
    
    prepare.output_fun <- function(spec.id,next_gen){
      
      spec.id$spec_event <- next_gen$n
      
      from.node <- spec.id$species_origin
      to.node <- spec.id$ID
      edge <- data.frame(cbind(from.node,to.node),spec_event = spec.id$spec_event,spp = spec.id$ID, alive = spec.id$alive)
      
      edge.length <- c()
      for (i in 1:(nrow(spec.id))){
        spec.origin <- spec.id$species_origin[i]
        edge.length[i] <- spec.id$generation[i]-spec.id$generation[which(spec.id$ID == spec.origin)]
      }
      
      stem.depth <- spec.id$generation
      
      edge$spp <- NA
      node <- 1
      spp <- 1
      
      a <- which(!is.na(all_species_matrix_f[,1]))
      b <- which(!is.na(all_species_matrix_f[,2]))
      a <- unique(c(a,b))
      edge$alive <- 0
      edge$alive[a] <- 1
      
      for (i in 1:nrow(edge)){
        
        if(edge$alive[i]==1){
          edge$node[i] = NA
          edge$spp[i] <- spp 
          spp <- spp + 1
          
        }else{
          edge$node[i] = node 
          node <- node + 1
        }
      }
      
      edge.f <- edge %>% 
        filter(spec_event > -999 | alive == 1) %>%
        mutate(ID = to.node)
      
      run<- unique(edge.f$from.node)
      
      n = -1
      for (i in run) {
        if (any(edge.f[, 'from.node'] == i)) {
          edge.f[which(edge.f[, 'from.node'] == i), 'from.node'] = n
          edge.f[which(edge.f[, 'to.node'] == i), 'to.node'] = n
          n = n - 1
        }
      }
      
      edge.f.f <- edge.f %>% 
        filter(!(to.node > 0 & alive == 0))
      
      
      filtered <- which(!(edge$to.node %in% edge.f.f$ID))
      
      if(identical(filtered,integer(0))){
        edge.out <- edge
        edge.length.out <- edge.length
        edge.length.out[1] <- 0
        stem.depth.out <- stem.depth
      }else{
        edge.out <- edge[-filtered,]
        edge.length.out <- edge.length[-filtered]
        edge.length.out[1] <- 0
        stem.depth.out <- stem.depth[-filtered]  
      }
      
      return(list(edge = edge.out,edge.length = edge.length.out,stem.depth = stem.depth.out))
    }
    
    make_phylo_fun <- function(output.list){
      edge.out <- output.list$edge
      edge.length.out <- output.list$edge.length
      stem.depth.out <- output.list$stem.depth
      
      edge.length.out[which(edge.out$alive == 1)] = sim_num.new - stem.depth.out[which(edge.out$alive == 1)]
      run<- unique(edge.out$from.node)
      edge.out.m = as.matrix(edge.out)
      n = -1
      
      for (i in run) {
        if (any(edge.out.m[, 'from.node'] == i)) {
          edge.out.m[which(edge.out.m[, 'from.node'] == i), 'from.node'] = n
          edge.out.m[which(edge.out.m[, 'to.node'] == i), 'to.node'] = n
          n = n - 1
        }
      }
      
      edge.out.m[which(edge.out.m[, c('from.node', 'to.node')] > 0)] = 1:length(edge.out.m[which(edge.out.m[, c('from.node', 'to.node')] > 0)])
      tip.label = edge.out.m[, 'to.node'][!is.na(edge.out[, 'spp'])]
      edge.only = edge.out.m[, c('from.node', 'to.node')]
      mode(edge.only) = "character"
      mode(tip.label) = "character"
      obj = list(edge = edge.only, edge.length = edge.length.out, tip.label = tip.label)
      class(obj) = "phylo"
      phylo.out = old2new.phylo(obj)
      return(phylo.out)
    }
    
    
    
    exp_seq.fun <- function(start,end,length){
      return(exp(seq(log(start), log(end), length.out = length)))
    }
    
    #probability of dispersal
    p_disp_trop <- sim_vec$disp_trop
    p_disp_temp <- sim_vec$disp_temp
    #probability of extinction
    p_xtc_trop <- sim_vec$xtc_trop
    p_xtc_temp <- sim_vec$xtc_temp
    #probability of speciation
    p_spec_trop <- sim_vec$spec_trop
    p_spec_temp <- sim_vec$spec_temp
    
    stop_max_spec <- stop_max_spec
    max_spec <- max_spec
    min_spec <- min_spec
    
    # run sim -----------------------------------------------------------------
    
    species_vector <- matrix(ncol = 2, nrow = 1)
    #colnames(species_vector) <- c("Trop","Temp")
    species_vector[1,] <- c(1,1)
    
    all_species_matrix <- matrix(NA, nrow = empty_row, ncol = 2)
    all_species_matrix[1,] <- species_vector
    
    #all_species_list <- list()
    #all_species_list[[1]] <- species_vector
    spec.id <- data.frame(species_origin = 0,
                          generation = 0,
                          source = "Trop",
                          ID = 1, 
                          alive = 1, 
                          speciated = NA
    )
    
    next_gen <- data.frame(n=-999,spec = 1)
    
    sim.num <- 1
    N.eqR <- N.eqR
    alive_species_count_trop <- numeric()  # For Trop species
    alive_species_count_temp <- numeric()  # For Temp species
    disp_fun_triggered <- FALSE
    
    repeat{
      
      n_spec <- nrow(spec.id)
      new.spec <- 1
      alive_sum <- sum(spec.id$alive)
      
      use_disp_fun <- alive_sum >= N.eqR
      
      if (alive_sum > stop_max_spec || alive_sum == 0 || sim.num > max.sim) {
        break
      }
      
      
      # if (alive_sum == 0 || sim.num > max.sim) {
      #      break
      #    } 
      #    
      for (j in 1:n_spec){
        
        species_vector <- all_species_matrix[j,]
        
        if (isTRUE(all(is.na(species_vector))) || sum(species_vector) == 0) {
          next 
        }
        
        if (use_disp_fun || disp_fun_triggered) {
          disp.ext_event <- disp_xtc.fun(species_vector)
          species_vector_new <- species_vector + disp.ext_event[1,] + disp.ext_event[2,]
          
          # Ensure this block always gets executed in future iterations
          disp_fun_triggered <- TRUE
        } 
        else {
          species_vector_new <- species_vector
          
        }
        
        if (sum(species_vector_new) ==0){
          all_species_matrix[j,] <- species_vector_new
          spec.id$alive[j] <- 0
          next
        }
        
        
        spec_event <- c(rbinom(1, 1, p_spec_trop),rbinom(1,1,p_spec_temp))
        speciation_event <- spec_event*species_vector_new
        
        if (any(speciation_event != 0)) {
          species_vector_new_matrix <- matrix(species_vector_new, nrow = 1, ncol = 2)
          
          all_species_matrix[j,] <- matrix(ncol = 2, nrow = 1, c(NA, NA))
          spec.id$alive[j] <- 0
          spec.id$speciated[j] = "Y"
          
          speciation_event_time <- speciation_event * round(runif(length(speciation_event)),5)
          
          new_species <- list()  # Temporarily store new species
          
          # This replaces multiple conditions checking speciation_event values
          new_species_counter <- 0
          
          if (sum(speciation_event) == 1){
            
            idx <- which(speciation_event == 1) %>% as.numeric()
            
            new_species_matrix <- matrix(ncol = 2, nrow = 1, data = 0)
            new_species_matrix[1, idx] <- 1
            
            #new_species <- append(new_species, list(species_vector_new_matrix, matrix(ncol = 2, nrow = 1, c(ifelse(idx==1,1,0), ifelse(idx==2,1,0)))))
            new_species <- append(new_species, list(species_vector_new_matrix, new_species_matrix))
            
            new_spec <- data.frame(species_origin = rep(j,2),
                                   generation = rep(sim.num+(speciation_event_time[idx]),2),
                                   source = ifelse(idx == 1, "Trop", "Temp"),
                                   ID = n_spec + new.spec + new_species_counter + (0:1),
                                   alive = rep(1,2),
                                   speciated = NA)
            
            spec.id <- bind_rows(spec.id,new_spec)
            
            
            next_gen[j,1] <- next_gen[j,1]+2
            
            new_gen <- data.frame(n = rep(-999,2), spec = n_spec + new.spec + new_species_counter + (0:1))
            
            next_gen <- bind_rows(next_gen,new_gen)
            
            new_species_counter <- new_species_counter + 2
            
          } else{
            
            sorted_times <- order(speciation_event_time)
            spec.id_prev <- spec.id[FALSE,]
            
            
            for (idx in sorted_times) {
              
              if (speciation_event[idx] != 0) {
                spec_origin_in_loop <- j
                
                
                #this is to remove the species that speciates twice
                if (new_species_counter>0){
                  
                  new_species[[new_species_counter-1]] <-  matrix(ncol = 2, nrow = 1, c(NA, NA))
                  spec.id_prev$alive[new_species_counter-1] <- 0
                  spec.id_prev$speciated[new_species_counter-1] <- "Y"
                  #spec_origin_in_loop <- j+(new_species_counter)
                  spec_origin_in_loop <- n_spec + new.spec + new_species_counter - 2
                  
                }
                
                new_species_matrix <- matrix(ncol = 2, nrow = 1, data = 0)
                new_species_matrix[1, idx] <- 1
                
                new_species <- append(new_species, list(species_vector_new_matrix,new_species_matrix))
                
                #new_species <- append(new_species, list(species_vector_new, matrix(ncol = 2, nrow = 1, c(ifelse(idx==1,1,0), ifelse(idx==2,1,0)))))
                
                new_spec <- data.frame(species_origin = rep(spec_origin_in_loop,2),
                                       generation = rep(sim.num+(speciation_event_time[idx]),2),
                                       source = ifelse(idx == 1, "Trop", "Temp"),
                                       ID = n_spec + new.spec + new_species_counter + (0:1),
                                       alive = rep(1,2),
                                       speciated = NA)
                
                
                
                spec.id_prev <- bind_rows(spec.id_prev,new_spec)
                
                next_gen[spec_origin_in_loop,1] <- next_gen[spec_origin_in_loop,1]+2
                
                new_gen <- data.frame(n = rep(-999,2), spec = n_spec + new.spec + new_species_counter + (0:1))
                
                next_gen <- bind_rows(next_gen,new_gen)
                
                new_species_counter <- new_species_counter + 2
                
              }
              
            }
            
            spec.id <- bind_rows(spec.id,spec.id_prev)
            
            
          }
          
          
          new.spec <- new.spec+ (sum(speciation_event) * 2)
          new_species <- do.call(rbind, lapply(new_species, function(x) as.matrix(x)))
          
          start <- nrow(spec.id)-(new_species_counter-1)
          end <- nrow(spec.id)
          
          all_species_matrix[start:end,] <- new_species 
          #all_species_list <- append(all_species_list,new_species)
          
          
        } else {
          all_species_matrix[j,] <- species_vector_new
          
        }  
        #now depending if there were or weren't speciation events:
        
        #no speciation events, keep the original matrix, and add a new generation 
        
        
      } 
      
      sim.num<- sim.num+1
      #alive_count_for_current_gen_trop <- sum(spec.id$alive == 1 & spec.id$source == "Trop")
      #alive_species_count_trop <- c(alive_species_count_trop, alive_count_for_current_gen_trop)
      
      #alive_count_for_current_gen_temp <- sum(spec.id$alive == 1 & spec.id$source == "Temp")
      #alive_species_count_temp <- c(alive_species_count_temp, alive_count_for_current_gen_temp)
      
    }
    
    
    sim_num.new <- sim.num
    
    if (nrow(spec.id)==1){
      
      #all_species_array =array(unlist(all_species_list),dim=c(nrow(all_species_list[[1]]),2,length(all_species_list)))
      #sum.arr <- apply(all_species_array, 1:2, sum,na.rm = T) %>% as.data.frame() %>% rename(Temperate =  V2, Tropic = V1) %>% reshape2::melt()
      #sum.arr$gen <- rep(1:(nrow(sum.arr)/2),2)
      
      return(list(out = "F", reason = "no species"))
      
      
    }else if (sum(spec.id$alive)<min_spec){
      
      #all_species_array =array(unlist(all_species_list),dim=c(nrow(all_species_list[[1]]),2,length(all_species_list)))
      #sum.arr <- apply(all_species_array, 1:2, sum,na.rm = T) %>% as.data.frame() %>% rename(Temperate =  V2, Tropic = V1) %>% reshape2::melt()
      #sum.arr$gen <- rep(1:(nrow(sum.arr)/2),2)
      
      
      return(list(out = "F",reason = "too few spec", spec.id = spec.id, sim.num = sim.num, total_spec = nrow(spec.id)))
      
    }else if (sum(spec.id$alive)>max_spec){
      
      return(list(out = "F",reason = "too many spec", spec.id = spec.id, sim.num = sim.num, total_spec = nrow(spec.id)))  
      
    }else{
      
      
      all_species_matrix_f <- all_species_matrix[1:nrow(spec.id),]
      
      
      
      
      # prepare output ----------------------------------------------------------
      
      
      output.list <- prepare.output_fun(spec.id,next_gen)
      
      # make phylo --------------------------------------------------------------
      
      
      
      phylo.sim <- make_phylo_fun(output.list)
      
      #if the phylogenetic tree doesn't work for some reason. 
      continue <- tryCatch({write.tree(phylo.sim)}, error = function(e){    
        
        continue <- "NO"
        return(continue)})
      
      if (continue == "NO"){
        
        return(list(out = "F", reason = "failed phylo"))
        
      }else{
        phylo.sim = read.tree(text = write.tree(phylo.sim))
        
        edge <- output.list$edge
        spp.alive <- data.frame(spp = edge$spp[(edge$alive==1)], ID = edge$to.node[edge$alive==1])
        
        alive_geo <- all_species_matrix_f[spp.alive$ID,]
        rownames(alive_geo) <- spp.alive$spp
        colnames(alive_geo) <- c("Trop","Temp") 
        
        phydist <- cophenetic.phylo(phylo.sim)
        mpd <-  ses.mpd(t(alive_geo),phydist, null.model = "taxa.labels" ,runs = 50, abundance.weighted = F)
        mpd.df <- c(MPD.z.1 = mpd$mpd.obs.z[1], MPD.z.2 = mpd$mpd.obs.z[2])
        
        return(list(out = "Y", sim_num = sim_num.new, spec.id = spec.id, next_gen = next_gen, all_species_array = all_species_matrix_f, mpd.df = mpd.df))
      }
    } 
  }
  
  
  sim_vec <- sim_mat[i,]
  
  #sim_mat <- data.frame(sim_mat_temp,sim_vec_trop)
  
  rep.cont <- 1
  max_runs <- 1  
  
  repeat {
    # Stop if we reach the maximum number of trials
    if (max_runs > max.reps) {
      break
    }
    
    if (rep.cont > reps) {
      
      break
    }
    
    if(run.bridge == T){
      sim.out <- sim_fun_bridge(sim_vec, stop_max_spec, min_spec) 
      
    } else {
      sim.out <- sim_fun_nobridge(sim_vec, stop_max_spec, min_spec) 
      
    }
    
    sim_out_list_par <- append(sim_out_list_par, list(sim.out))
    
    if (sim.out$out == "F") {
    } else {
      
      rep.cont <- rep.cont + 1
    }
    
    max_runs <- max_runs + 1
  }
  
  
  sim.vec.list <- list(param = sim_vec)

  col.choose = ifelse(run.bridge == T, 3,2)
  
  mpd.mat <- matrix(nrow = reps, ncol = col.choose)
  sim.c <- as.numeric(sim_vec)
  for (k in 1:reps){
    mpd.rep <- sim_out_list_par[[k]]$mpd.df
    if (is.null(mpd.rep)){
      mpd.rep <- rep(NA,col.choose)
    }
    mpd.mat[k,1:col.choose] <- mpd.rep
  }
  
  mpd.df <- mpd.mat %>% as_data_frame()
  
  output.df <- mpd.df %>%
    summarise_all(list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)))
  output.df$MPDdif <- output.df$V2_mean - output.df$V1_mean
  #colnames(output.df)[1:6] <- c("MPD1", "MPD2","MPD3","sd1","sd2","sd3")
  output.df <- cbind(output.df,sim_vec, row.names = NULL)
  
  mpd.list <- list(mpd.df = mpd.df,output = output.df)
  sim_out_list_par_out <- c(list(runs = sim_out_list_par), sim.vec.list, mpd.list)
  return(sim_out_list_par_out)
  
}
