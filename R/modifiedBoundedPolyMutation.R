myBoundedPolyMutation <-
  function(parent_chromosome,lowerBounds,upperBounds,mprob,mum){
    popSize=nrow(parent_chromosome);
    varNo=ncol(parent_chromosome);
    child <- parent_chromosome;
    for (i in 1:popSize) {
      for (j in 1:varNo){
        y = child[i,j];
        # if the random probability is less than mprob, then mutate that variable
        if (runif(1) < mprob) {
          yl = lowerBounds[j];
          yu = upperBounds[j];
          if (y > yl) {
            rnd = runif(1);
            mut_pow = 1.0/(mum + 1.0);
            if (rnd < 0.5){
              delta = (y-yl);
              val = 2.0*rnd
              deltaq = val^mut_pow - 1.0;
              y = y + deltaq*delta
            } else {
              delta = (yu-y);
              val = 2.0*(1.0-rnd)
              deltaq = 1.0 - val^mut_pow;
              y = y + deltaq*delta
            }

            if (y > yu) {
              y = yu;
            } else if (y < yl) {
              y = yl;
            }
            child[i,j] = y;
          } else { # y <= yl
            xy = runif(1);
            child[i,j] = yl + xy*(yu-yl);
          }
        } # runif(1) > mprob, do not perform mutation
      } # next j
    } # next i
    return(child);
  }
