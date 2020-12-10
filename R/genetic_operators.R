uniformXover <- function(parent_chromosome,
                         lowerBounds=NULL,
                         upperBounds=NULL,
                         cprob = 0.9){
  chromosome_length <- ncol(parent_chromosome)

  nPair <- floor(nrow(parent_chromosome)/2)
  if(is.null(lowerBounds)){
    lowerBounds <- rep(0,chromosome_length)
  }
  if(is.null(upperBounds)){
    upperBounds <- rep(1,chromosome_length)
  }
  offspring <- NULL
  altoffspring_collection <- NULL
  for(parentPair in 1:nPair){
    thisoffspring <- parent_chromosome[parentPair*2-1,]
    altoffspring <- parent_chromosome[parentPair*2,]
    for(i in 1:chromosome_length){
      if(runif(1)<cprob){
        thisoffspring[i] <- parent_chromosome[parentPair*2,i]
        altoffspring[i] <- parent_chromosome[parentPair*2-1,i]
      }
    }
    offspring <- rbind(offspring,matrix(thisoffspring,nrow=1))
    altoffspring_collection <- rbind(altoffspring_collection,matrix(altoffspring,nrow=1))
  }
  offspring <- rbind(offspring,altoffspring_collection)

  return(offspring)
}


truncnormMutation <- function(parent_chromosome,
                              lowerBounds=NULL,
                              upperBounds=NULL,
                              mprob = NULL,
                              sigma = NULL){
  chromosome_length <- ncol(parent_chromosome)
  nOffspring <- nrow(parent_chromosome)
  if(is.null(mprob)){
    mprob <- 1/chromosome_length
  }
  if(is.null(lowerBounds)){
    lowerBounds <- rep(0,chromosome_length)
  }
  if(is.null(upperBounds)){
    upperBounds <- rep(1,chromosome_length)
  }
  if(is.null(sigma)){
    sigma <- rep(1,chromosome_length)
  }
  var_range <- upperBounds-lowerBounds
  offspring <- NULL
  for(offspringIndex in 1:nOffspring){
    thisoffspring <- parent_chromosome[offspringIndex,]
    for(i in 1:chromosome_length){
      if(runif(1)<mprob){
        thisoffspring[i] <- truncnorm::rtruncnorm(1,
                                                  a=lowerBounds[i],
                                                  b=upperBounds[i],
                                                  mean=thisoffspring[i],
                                                  sd=sigma)
      }
    }
    offspring <- rbind(offspring,matrix(thisoffspring,nrow=1))
  }
  return(offspring)
}

# if mirrorred, 2 offspring are created
orthogonal_sampling_mutation <- function(parent_chromosome,
                                         lowerBounds=NULL,
                                         upperBounds=NULL,
                                         mprob = NULL,
                                         sigma = NULL,
                                         nOffspring=1,
                                         nDirection=1,
                                         searchDir=NULL){
  reference <- parent_chromosome
  chromosome_length <- ncol(reference)
  if(is.null(lowerBounds)){
    lowerBounds <- rep(0,chromosome_length)
  }
  if(is.null(upperBounds)){
    upperBounds <- rep(1,chromosome_length)
  }
  if(is.null(sigma)){
    sigma <- rep(1,nDirection)
  }
  var_range <- upperBounds-lowerBounds
  offspring <- NULL

  # col major
  if(is.null(searchDir)){
    o <- pracma::gramSchmidt(matrix(runif(chromosome_length*nDirection),
                                    nrow=chromosome_length))$Q#[,1:nDirection]
  }else{
    o <- pracma::gramSchmidt(searchDir)$Q[,1:nDirection]
  }
  # orthogonal_direction <- abs(orthogonal_direction)

  for(offspringIndex in 1:nOffspring){
    thisoffspring <- reference[offspringIndex,]
    # check the bound w.r.t. the orthogonal direction
    for(dirIndex in 1:nDirection){
      dist_to_ubound <- upperBounds-thisoffspring
      dist_to_lbound <- lowerBounds-thisoffspring

      transformed_dist_to_lbound <- (dist_to_lbound/o[,dirIndex])
      transformed_dist_to_ubound <- (dist_to_ubound/o[,dirIndex])
      transformed_dist <- c(transformed_dist_to_lbound,transformed_dist_to_ubound)
      positives <- which(transformed_dist>=0)
      negatives <- which(transformed_dist<0)

      allowed_range <- max(transformed_dist[negatives])
      allowed_range <- c(allowed_range,min(transformed_dist[positives]))

      ortho_range <- allowed_range[2]-allowed_range[1] # multiplier range

      # testUbound <- append(testUbound,shortest_dist_to_ubound)
      # testLbound <- append(testLbound,shortest_dist_to_lbound)
      # truncated normal distribution is used to prevent
      mutation_scale <- truncnorm::rtruncnorm(1,
                                              a=allowed_range[1],
                                              b=allowed_range[2],
                                              mean=0,
                                              sd=sigma[dirIndex])
      thisoffspring <- thisoffspring+ o[,dirIndex]*mutation_scale
    }

    # print(thisoffspring+orthogonal_direction%*%testUbound)

    # thisoffspring <- thisoffspring*var_range
    offspring <- rbind(offspring,matrix(thisoffspring,nrow=1))
  }
  # if(mirrorred){
  offspring_dist <- offspring - reference

  mirrorred_offspring <- reference - offspring_dist
  if(all(mirrorred_offspring<upperBounds) && all(mirrorred_offspring>lowerBounds)){
    mirrorred_offspring <- mirrorred_offspring
  }else{
    mirrorred_offspring <- NULL
  }
  # }
  return(list(offspring=offspring,mirror=mirrorred_offspring,searchDir=o))
}
