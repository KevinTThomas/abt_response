mem.profile.linux <- function() {
  # Split system information
  info_split <- strsplit(system("free", intern = TRUE), split = " ")
  
  # Remove "Mem:" and "Swap:"
  info_split <- lapply(info_split, function(x){gsub("Mem:", "", x)})
  info_split <- lapply(info_split, function(x){gsub("Swap:", "", x)})
  
  # Get actual values
  info_split <- lapply(info_split, function(x){x[x != ""]})
  
  # Bind values
  info_split <- do.call(rbind, info_split[1:2])
  df <- as.data.frame(info_split[2,,drop=FALSE])
  names(df) <- info_split[1,]
  df <- as.data.frame(lapply(df, as.double))
  
  return(df)
}
