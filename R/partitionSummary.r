tosummaryMat = function(x){
	#Turns the return from partionSummary into a matrix
	return(data.frame(id =x$levels_out,
		n = x$n_vec))
}
partitionSummary = function(outcome, part1, part2 = NULL){
  ##outcome should be numeric
  #part1, part2, are partitioning factors
  #returns:  levels, labels, n_vec, mean_vec
  #This function doesnt handle names that well,
  # and the (few) ways we partition the data must be mmanually labeled
  #Can optionally add dg pcf and group vecs to obtain Draper-Gittoes scores
  #Optimized by Lucy Greenberg
if (is.null(part1)){return(NULL)}
  #Both vectors need to be partitions
 part1 = as.factor(part1)
 if (!is.null(part2)){
   stopifnot(length(part1) == length(part2))
   part2 = as.factor(part2)
   }
   
  #If two partition vectors specified, run through data
  # and list all combinations of levels (e.g. inst-subst)
  if (is.null(part2)){
    part_mat = data.frame(levels(part1))
  } else{
    part_mat = data.frame(part1 = NULL, part2 = NULL)
    for (p1 in levels(part1)){
      p2_vec = sort(unique(part2[part1==p1]))
      part_mat = rbind(part_mat,
                       cbind( as.character(p1), as.character(p2_vec)))
    }
    #names(part_mat) = c(names(part1), names(part2))
  }
  
  #vector partition levels of length q
  if (is.null(part2)){
    levels_out = levels(part1)
  } else{
    levels_out = apply(part_mat, 1, paste, collapse = '-')
  }

  #Vector partitions of length n
  if (is.null(part2)){
    part_vec = part1
  } else {
    part_vec = factor(apply(cbind(as.character(part1),as.character(part2)), 1, paste, collapse = '-'),
		levels = levels_out) #Need to manually set levels or the order changes
  }
  v_mat = cbind(as.character(part1),as.character(part2))
 
  #Extract variables
  q = length(levels_out)
  n = length(outcome)
  stopifnot( n == length(part_vec)) #Partion vector same length as outcome

  #Partition indicator matrix (n x q)
  #Useful for a lot of linear algebra shortcuts
  #ex: o_mean =  t(ind_mat) %*% ind_mat / n_vec
  ind_mat = Matrix::Matrix(0, nrow = n, ncol = q)

  ind_mat[cbind(1:n, as.numeric(part_vec))] = 1

  #Statistics for each partition
  # n_vec = sapply(levels_out, function(li) sum(part_vec == li))
  # o_mean = sapply(levels_out, function(li) mean(outcome[part_vec == li]))
  temp <- data.frame(table(part_vec))
      n_vec<-temp[,2]
      names(n_vec)<-temp[,1]
      
      temp2<-aggregate(outcome,by=list(part_vec),mean)
      o_mean<-temp2[,2]
      names(o_mean)<-temp2[,1]

  return(list(levels_out = levels_out,
              n_vec = c(n_vec),
              o_mean = c(o_mean),
              ind_mat = ind_mat,
              part_mat = part_mat
			  ))
}
