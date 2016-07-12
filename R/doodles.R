`terminal.twigs` <-
  function(topology){
    n<-length(topology);
    return(seq(max(topology)+1, n));
  }


`parse.tree` <-
  function(topology, times){
    term<-terminal.twigs(topology);
    N <- length(term);
    dm<-distance.matrix(topology, times);
    bt<-branch.times(topology, times);
    pt<-list(N=N, term=term, dm=dm, bt=bt);
    return(pt);
  }


`distance.matrix` <-
  function (topology, times) {
    term <- terminal.twigs(topology);
    N <- length(term);
    dm <- matrix(data=0,nrow=N,ncol=N);
    dm[1,1] <- 0;
    for (i in 2:N) {
      pedi <- pedigree(topology,term[i]);
      for (j in 1:(i-1)) {
        pedj <- pedigree(topology,term[j]);
        for (k in 1:length(pedi)) {
          if (any(pedj == pedi[k])) break;
        }
        dm[j,i] <- dm[i,j] <- (times[term[i]]-times[pedi[k]]) + (times[term[j]]-times[pedi[k]]);
      }
      dm[i,i] <- 0;
    }
    return(dm);
  }

`pedigree` <-
  function (topology, k) {
    p <- k;
    k <- topology[k];
    while (k != 0) {
      p <- c(p, k);
      k <- topology[k];
    }
    return(p);
  }


`branch.times` <-
  function (topology, times) {
    term <- terminal.twigs(topology);
    N <- length(term);
    bt <- matrix(data=0,nrow=N,ncol=N);
    bt[1,1] <- times[term[1]];
    for (i in 2:N) {
      pedi <- pedigree(topology,term[i]);
      for (j in 1:(i-1)) {
        pedj <- pedigree(topology,term[j]);
        for (k in 1:length(pedi)) {
          if (any(pedj == pedi[k])) break;
        }
        bt[j,i] <- bt[i,j] <- times[pedi[k]];
      }
      bt[i,i] <- times[term[i]];
    }
    return(bt);
  }


N = 86
num.prob<-matrix(data=0, nrow=N, ncol=N)


# pt<-parse.tree(topology, times);
pt <- parse.tree(input$ancestor, input$time)
ta<-pt$bt

alpha.est <- 1

for(p in 1:N)
{
  for(q in 1:N)
  {
    if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
  }
}
num.prob<-matrix(data=1, nrow=N, ncol=N)
microbenchmark(
  a = test2(ta, alpha.est = 1),
  ac = test2c(ta, alpha.est = 1),
  applyfunc = (num.prob.1 <- matrix(vapply(ta.vec, num.prob.fc, FUN.VALUE = 0), ncol=N, nrow=N))
)

test2 <- function(ta, num.prob, alpha.est){
  num.prob <- matrix(1, nrow=86, ncol=86)
  for(p in 1:N) {
    for(q in 1:N) {
      if(ta[q,p]==0) num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
    }
  }
return(num.prob)
}
test2c <- cmpfun(test2)

ta.vec <- as.vector(ta)
num.prob.f <- function (x){ if (x != 0) return ((1-exp(-alpha.est*x))/(alpha.est*x)) else return(1) }
num.prob.fc <- cmpfun(num.prob.f)
num.prob.1 <- matrix(vapply(ta.vec, num.prob.f, FUN.VALUE = 0), ncol=N, nrow=N)


identical (num.prob.1, num.prob)

as.matrix()

aljaball()


