#adopted from Jon Lefcheck 'jslefche' - https://gist.github.com/jslefche/47a48fad73d4d5a82005

# Function estimateR() takes a species (rows)-by-community (columns) abundance matrix, 
# order q, specifying species richness (q = 0, default), Shannon entropy (q = 1), inverse Simpson index (q = 2)
# a minimum level of min.coverage for estimation (if specified), 
# a minimum level of min.sample for rarefaction (if rare == TRUE),
# whether a histogram of the coverage values should be plotted (plotit = TRUE),
# whether a progress bar should be shown (.progressBar = TRUE)
# and additional arguments to iNEXT() (...)

# The function returns a data.frame containing observed richness, coverage at the observed level of richness,
# estimated richness, Chao1 index (if chao1 = TRUE) and rarefied richness (individual-based, if rare = TRUE)

estimateR = function(
  comm, min.coverage = NULL, min.sample = NULL, chao1 = TRUE, rare = TRUE, 
  plotit = TRUE, .progressBar = TRUE, ...) {
  
  # If comm is not a list, turn each row into an entry in a list
  if(!is.list(comm)) comm = as.list(data.frame(comm))
  
  # Set progress bar
  if(.progressBar == TRUE)  pb = txtProgressBar(min = 0, max = length(comm), style = 3) else pb = NULL
  
  # Run iNEXT function to estimate sample coverage 
  iNEXT.list = lapply(1:length(comm), function(i) {
    
    if(sum(comm[[i]]>0) > 1) {
      
      try.error = try(iNEXT(as.numeric(comm[[i]])))
      
      if(!is.null(pb)) setTxtProgressBar(pb, i)
      
      if(class(try.error) == "try-error") NA else try.error
      
    } else NA
    
  } )
  
  # Close progress bar
  if(!is.null(pb)) close(pb)  
  
  # Isolate errors and attempt to fix
  if(any(unlist(lapply(iNEXT.list, is.na)))) {
    
    replace.list = lapply(which(lapply(iNEXT.list, function(i) all(is.na(i))) == TRUE), function(i) { 
      
      replace.list = list(
        
        data.frame(n = sum(comm[[i]]), S.obs = sum(comm[[i]] > 0), C.hat = NA, 
                   f1 = NA, f2 = NA, f3 = NA, f4 = NA, f5 = NA,
                   f6 = NA, f7 = NA, f8 = NA, f9 = NA, f10 = NA),
        
        matrix(c(rep(sum(comm[[i]] > 0), 6), rep(0, 9)), nrow = 3, ncol = 5, 
               dimnames = 
                 list(c("Species Richness", "Expontential Entropy", "Inverse Simpson"), 
                      c("Observed", "Estimator", "Est_s.e.", "95% Lower", "95% Upper"))),
        
        data.frame(m = sum(comm[[i]]), method = "observed", order = 0, 
                   qD = ifelse(sum(comm[[i]]) == 0, NA, sum(comm[[i]] > 0)),
                   qD.95.LCL = 0, qD.95.UCL = 0, SC = NA, SC.95.LCL = NA, SC.95.UCL = NA)
        
      )
      
      names(replace.list) = c("DataInfo", "BasicIndex", "Accumulation")
      
      return(replace.list)
      
    } )
    
    # Replace in master list
    # Index appropriate entry in master list
    to.replace = which(lapply(iNEXT.list, function(i) all(is.na(i))) == TRUE)
    
    for(i in 1:length(replace.list)) {
      
      iNEXT.list[[to.replace[i]]] = replace.list[[i]]
      
    }
    
  }
  
  # Convert list of lists into a single list
  iNEXT.list = list(
    
    # $DataInfo
    do.call(rbind, lapply(1:length(iNEXT.list), function(i) return(iNEXT.list[[i]]$DataInfo) )),
    
    # $BasicIndex
    do.call(rbind, lapply(1:length(iNEXT.list), function(i) return(iNEXT.list[[i]]$BasicIndex) )),
    
    # $Accumulation
    do.call(list, lapply(1:length(iNEXT.list), function(i) return(iNEXT.list[[i]]$Accumulation) ))
    
  )
  
  names(iNEXT.list) = c("DataInfo", "BasicIndex", "Accumulation")
  
  # Get minimum level of coverage across all samples (unless user-supplied)
  if(is.null(min.coverage)) {
    
    all.coverage = sapply(iNEXT.list$Accumulation, function(x) x[x$method == "observed", "SC"])
    
    if(plotit == TRUE) {
      
      if(is.null(min.sample)) par(mfrow = c(1, 2)) 
      
      hist(all.coverage[!is.na(all.coverage)], xlab = "Sample coverage", main = "")
      
    }
    
    min.coverage = min(sapply(iNEXT.list$Accumulation, function(x)
      
      x[x$method == "observed", "SC"]), na.rm = TRUE)
    
    print(paste("The minimum level of coverage across all samples was ", 100 * min.coverage, "%", sep = ""))
    
  }
  
  # Estimate richness at the minimum level of coverage
  estimate.S = sapply(iNEXT.list$Accumulation, function(x)
    
    # If SC is NA, return S.obs
    if(all(is.na(x$SC))) x$qD else
      
      # If sample coverage is complete (i.e., 1) for all levels of m, reported observed richness
      if(all(x$SC == 1)) x[x$method == "observed", "qD"] else
        
        # If user-supplied coverage is less than maximum coverage sampled, return NA
        if(min.coverage > max(x$SC)) NA else 
          
          # Else return linear approximation
          approx(x = x$SC, y = x$qD, xout = seq(0, 1, 0.001), rule = 2)$y[min.coverage * 1000]
    
  )
  
  # Return results in a data.frame
  diversity.df = data.frame(
    observed.S = iNEXT.list$DataInfo[, 2],
    observed.Chat = iNEXT.list$DataInfo$C.hat,
    observed.N = iNEXT.list$DataInfo$n,
    estimated.S = estimate.S
  )
  
  # Get Chao1 estimate if Chao1 == TRUE
  if(chao1 == TRUE) 
    
    diversity.df = data.frame(
      diversity.df, 
      Chao1 = matrix(iNEXT.list$BasicIndex[, 2], nrow = 3)[1, ]
    )
  
  # Get rarified richness if rare == TRUE
  if(rare == TRUE) {
    
    if(is.null(min.sample)) {
      
      # Calculate minimum number of individuals observed in a given sample
      all.m = sapply(iNEXT.list$Accumulation, function(x) 
        
        if(all(is.na(x$SC))) NA else x[x$method == "observed", "m"]
        
      )
      
      if(plotit == TRUE) hist(all.m[!is.na(all.m)], xlab = "Number of Individuals", main = "")
      
      min.sample = min(all.m, na.rm = TRUE)
      
      print(paste("The minimum number of individuals across all samples was ", min.sample, sep = ""))
      
    }
    
    # Interpolate backwards to estimate richnesss for that number of individuals across all samples
    rare.S = sapply(iNEXT.list$Accumulation, function(x)
      
      if(nrow(x) == 1) x$qD else
        
        if(min.sample > max(x$m)) NA else
          
          approx(x = x$m, y = x$qD, xout = seq(0, max(x$m), 1), rule = 2)$y[min.sample]
      
    )
    
    diversity.df = data.frame(
      diversity.df,
      rarefied.S = rare.S
    )
    
  }
  
  # Return data.frame
  return(diversity.df)
  
}