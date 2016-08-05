## Extra functions

setClass("slouch", contains = "list")

fitch<-function(tree.data, niche, deltran=FALSE, acctran=FALSE, root=NULL){
  niche<-as.factor(niche)
  anc.char<-as.character(tree.data$ancestor)
  anc.num<-as.numeric(anc.char)
  anc.num[1]<-0
  node.char<-as.character(tree.data$nodes)
  node.num<-as.numeric(node.char)
  niche.num<-as.numeric(niche)
  niches<-levels(niche)
  num.code<-1:length(levels(niche))
  translator<-cbind(niches, num.code)
  obs.char<-unique(niche.num[!is.na(niche.num)])
  N<-length(node.num)


  tree.matrix<-matrix(data=c(anc.num, node.num, niche.num), ncol=3)
  colnames(tree.matrix)<-c("Ancestors", "Nodes", "States")
  Nint<-max(tree.matrix[,1])
  Nnodes<-max(tree.matrix[,2])



  downpass.node.states<-list()
  cost<-0

  # place tip character states in downpass list

  for(i in 1:Nnodes)
  {
    downpass.node.states[[i]]<-tree.matrix[i,3]
  }


  check.node.order<-tree.matrix[,2]-tree.matrix[,1]
  if(any(check.node.order<0))
    traversal<-2:Nint else traversal<-Nint:1

  # Fitch dwonpass algorithm for internal nodes

  for(i in traversal)
  {

    children<-which(tree.matrix[i,2]==tree.matrix[,1])


    if(length(children)==2)
    {
      child1.state<-downpass.node.states[[children[1]]];
      child2.state<-downpass.node.states[[children[2]]];
      downpass.node.states[[i]]<-intersect(child1.state, child2.state); 	 }

    # need to extend this for cases greater than tritomies;

    if(length(children) > 2)
    {
      child1.state<-downpass.node.states[[children[1]]];
      child2.state<-downpass.node.states[[children[2]]]; 	  child3.state<-downpass.node.states[[children[3]]];
      tmp<-intersect(child1.state, child2.state)
      downpass.node.states[[i]]<-intersect(tmp, child3.state); 	 }

    # if no intersect of states, use union and update costs for state change;

    if(length(downpass.node.states[[i]])==0)
    {
      if(length(children)==2)
      {
        downpass.node.states[[i]]<-union(child1.state, child2.state);
        cost = cost+1;
      }
      if(length(children)>2)
      {
        tmp<-union(child1.state, child2.state);
        downpass.node.states[[i]]<-union(tmp, child3.state);
        cost = cost+1;
      }
    }
  }

  # final traversal for root node if internal node order is reversed

  if(any(check.node.order<0))
  {
    children<-which(tree.matrix[1,2]==tree.matrix[,1])

    if(length(children)==2)
    {
      child1.state<-downpass.node.states[[children[1]]];
      child2.state<-downpass.node.states[[children[2]]];
      downpass.node.states[[1]]<-intersect(child1.state, child2.state); 	 }

    # need to extend this for any number of polytomies;

    if(length(children) > 2)
    {
      child1.state<-downpass.node.states[[children[1]]];
      child2.state<-downpass.node.states[[children[2]]]; 	  child3.state<-downpass.node.states[[children[3]]];
      downpass.node.states[[1]]<-intersect(c(child1.state, child2.state), child3.state); 	 }

    if(length(downpass.node.states[[1]])==0)
    {
      if(length(children)==2)
      {
        downpass.node.states[[1]]<-union(child1.state, child2.state);
        cost = cost+1;
      }
      if(length(children)>2)
      {
        tmp<-union(child1.state, child.state)
        downpass.node.states[[1]]<-union(tmp, child3.state);
        cost = cost+1;
      }
    }
  }  # end of reverse tree traversal with root done separately;


  # Fitch up pass algorithm
  # the downpass is not guaranteed to give the optimal
  # node states (as defined by minimum number of changes)
  # for this we need the slighly more complex Fitch up-pass
  # algorithm;

  if(any(check.node.order<0))
    traversal<-2:Nint else traversal<-Nint:2
  pre.traversal<-rev(traversal)


  # start with assuming that final state set for the root
  # is the same as the downpass state set for the root
  # note that if the root is ambiguous, we need to set it to a single state for the algorithm to work

  finalpass.node.states<-downpass.node.states



  if(!is.null(root)) finalpass.node.states[[1]]<-root

  if(length(finalpass.node.states[[1]]) >=2)
  {
    message("There is an ambiguity at the root, as given below")
    print(as.numeric(finalpass.node.states[[1]]))
    message("One needs to set this to one of the states using root = state in the function call before attempting the deltran or acctran reconstructions")
  }

  for(i in pre.traversal)
  {

    #parent<-which(tree.matrix[i,1]==tree.matrix[,2]);]
    parent<-tree.matrix[i,2]
    ancestor<-tree.matrix[i,1]
    children<-which(tree.matrix[,1]==tree.matrix[i,2]);

    # set final pass state for node i to the intersect of node i with its immediate ancestor's state

    finalpass.node.states[[parent]]<-intersect(downpass.node.states[[parent]], finalpass.node.states[[ancestor]]);

    if(setequal(finalpass.node.states[[parent]], finalpass.node.states[[ancestor]])==F)
    {
      # conditions for bifurcating

      if(length(children)==2)
      {
        child1.state<-downpass.node.states[[children[1]]];
        child2.state<-downpass.node.states[[children[2]]];
        if(length(intersect(child1.state, child2.state))!=0)
        {
          finalpass.node.states[[parent]]<-union(downpass.node.states[[parent]], (intersect(finalpass.node.states[[ancestor]], union(child1.state, child2.state))));
        }
        if(length(intersect(child1.state, child2.state))==0)
        {
          finalpass.node.states[[i]]<-union(downpass.node.states[[parent]], finalpass.node.states[[ancestor]])
        }

        if(deltran==TRUE){
          if(length(finalpass.node.states[[i]]) >=2)  finalpass.node.states[[i]]= finalpass.node.states[[ancestor]]}

        if(acctran==TRUE){
          if(length(finalpass.node.states[[i]]) >=2)  finalpass.node.states[[i]]= setdiff(finalpass.node.states[[parent]],finalpass.node.states[[ancestor]]) }


      }       # end of bifurcaring condition


      # condition for polytomy

      if(length(children)>2)
      {
        child1.state<-downpass.node.states[[children[1]]];
        child2.state<-downpass.node.states[[children[2]]];
        child3.state<-downpass.node.states[[children[3]]]

        if(length(intersect(child3.state, intersect(child1.state, child2.state)))!=0)
        {
          finalpass.node.states[[parent]]<-union(downpass.node.states[[parent]], (intersect(finalpass.node.states[[ancestor]], union(union(child3.state,  child1.state), child2.state))));
        }
        if(length(intersect(child3.state, intersect(child1.state, child2.state)))==0)
        {
          finalpass.node.states[[i]]<-union(downpass.node.states[[parent]], finalpass.node.states[[ancestor]])
        }

      }     # end of polytomy condition

    } # end of final parent state not equal ancestor state


  } # end of post tree traversal

  # enumerate all posible optimal state
  # reconstructions in dataframe (ouchtree) format
  # and reconstitue niches as characters;

  return(list(treematrix=tree.matrix, Final.states=finalpass.node.states, downpass.cost=cost, niche.code=translator))
}

make.states<-function(pars.object){
  x<-pars.object$Final.states
  n<-length(x)
  n.states<-length(pars.object$niche.code[,1])

  #which nodes are ambiguous

  count<-0

  ambig<-NA
  for(i in 1:n)
  {
    if(length(x[[i]])>=2)
    {
      count=count+1;
      ambig<-c(ambig,i)
    }
  }
  ambig=ambig[-1]
  n.ambig<-length(ambig)

  # choose first character state

  tmp<-matrix(data=NA, nrow=n, ncol=1)
  for(i in 1:n)
  {
    tmp[i,1]<-x[[i]][1]
    for(j in 1:n.states)
    {
      if (tmp[i,1]==j) tmp[i,1]<-pars.object$niche.code[,1][[j]]
    }
  }

  # encode ambiguous states as character ambiguous

  if(n.ambig!=0)
  {
    for(i in 1:n)
    {
      for(j in 1:n.ambig)
      {
        if(i==ambig[j]) tmp[i,1]="ambiguous"
      }
    }
  }

  return(as.factor(tmp))
}

n.ambig<-function(pars.states)
{
  states<-pars.states$Final.states
  N<-length(states)
  count<-0
  ambig<-NA
  for(i in 1:N)
  {
    if(length(states[[i]]) >=2)
    {
      count=count+1;
      ambig<-c(ambig, i)
    }

  }
  return(list(N.ambiguous=count, ambiguous.nodes=ambig[-1]))
}


`arrange.tree` <-
  function (root, topology) {
    k <- which(topology==root);
    n <- length(k);
    reltree <- rep(0,length(topology));
    reltree[root] <- 0.5;
    p <- list(NULL);
    if (n > 0) {
      m <- rep(0,n);
      for (j in 1:n) {
        p[[j]] <- arrange.tree(k[j],topology);
        m[j] <- length(which(p[[j]] != 0));
      }
      cm <- c(0,cumsum(m));
      for (j in 1:n) {
        reltree <- reltree + (cm[j]/sum(m))*(p[[j]] != 0) + (m[j]/sum(m))*p[[j]];
      }
    }
    return(reltree);
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


`epochs` <-
  function (topology, times, term) {
    N <- length(term);
    e <- vector(length=N,mode="list");
    for (k in 1:N) {
      p <- pedigree(topology,term[k]);
      e[[k]] <- times[p];
    }
    return(e);
  }


`make.tree` <-
  function(n.tips, stretch=NULL){
    n.nodes<-n.tips*2-1
    int.nodes<-n.nodes-n.tips
    species<-c(rep(NA, times=int.nodes), letters[1:n.tips])
    generic<-c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 48, 49, 49, 50, 50, 51, 51, 52, 52, 53, 53, 54, 54, 55, 55, 56, 56, 57, 57, 58, 58, 59, 59, 60, 60, 61, 61, 62, 62, 63, 63)
    if(n.tips==4)
    {
      ancestor<-generic[1:n.nodes]
      time<-c(0, 0.5, 0.5, 1, 1, 1, 1)
    }
    if(n.tips==6)
    {
      ancestor<-c(0, 1, 2, 2, 1, 3, 3, 4, 4, 5, 5)
      time<-c(0, 1/3, 2/3, 2/3, 2/3, rep(1, times=n.tips))
    }
    if(n.tips==8)
    {
      ancestor<-generic[1:n.nodes]
      time<-c(0, 1/3, 1/3, 2/3, 2/3, 2/3, 2/3, rep(1, times=n.tips))
    }
    if(n.tips==10)
    {
      ancestor<-c(0, 1, 1, 2, 3, 3, 4, 4, 2, generic[9:n.nodes])
      time<-c(0, 0.25, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75, rep(1, times=n.tips) )
    }
    if(n.tips==12)
    {
      ancestor<-c(0, 1, 2, 2, 1, generic[6:n.nodes])
      time<-c(0, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, rep(1, times=n.tips))
    }
    if(n.tips==14)
    {
      ancestor<-c(0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 3, generic[14:n.nodes])
      time<-c(0, 0.25, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, rep(1, times=n.tips))}
    if(n.tips==16)
    {
      ancestor<-generic[1:n.nodes]
      time<-c(0, 0.25, 0.25, rep(0.5, times=4), rep(0.75, times=8), rep(1, times=n.tips))
    }
    if(n.tips==18)
    {
      ancestor<-c(0, 1, 1, 3, 2, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 4, generic[18:n.nodes])
      time<-c(0, 0.2, 0.2, 0.4, rep(0.6, times=4), rep(0.8, times=9), rep(1, times=n.tips))
    }
    if(n.tips==20)
    {
      ancestor<-c(0, 1, 1, 2, 3, 3, 4, 4, 2, generic[10:n.nodes])
      time<-c(0, 0.2, 0.4, 0.4, rep(0.6, times=5), rep(0.8, times=10), rep(1, times=n.tips))
    }
    if(n.tips==32)
    {
      ancestor<-	generic[1:n.nodes]
      time<-c(0, 0.2, 0.2, rep(0.4, times=4), rep(0.6, times=8), rep(0.8, times=16), rep(1, times=n.tips))
      species<-c(rep(NA, times=int.nodes), 1:32)
    }
    if(n.tips==64)
    {
      ancestor<-generic[1:n.nodes]
      time<-c(0, 1/6, 1/6, rep(2/6, times=4), rep(0.5, times=8), rep(4/6, times=16), rep(5/6, times=32), rep(1, times=n.tips))
      species<-c(rep(NA, times=int.nodes), 1:64)
    }


    if(!is.null(stretch)) time<-time^stretch
    tree<-data.frame(ancestor, time, species)
    return(tree)
  }


`no.me.sigma.X.estimate` <-
  function (predictor, topology, times) {
    pt <- parse.tree(topology,times);
    n <- pt$N;
    v <- pt$bt;
    w <- matrix(data=1,nrow=pt$N,ncol=1);
    dat <- predictor[!is.na(predictor)];
    beta<-solve(t(w)%*%solve(v)%*%w)%*%(t(w)%*%solve(v)%*%dat)
    e<-dat-beta
    theta <- beta
    sigma <- sqrt((e %*% solve(v,e))/(n-1));
    dim(sigma) <- 1;
    return(list(as.numeric(theta), as.numeric(sigma)));
  }

`np.regression` <-
  function(response, me.response, predictor, me.predictor, convergence=NULL){
    if(is.null(convergence)) convergence=0.000001
    Y <- response[!is.na(response)];
    N <- length(Y);
    pred<-data.frame(predictor)
    me.pred<-data.frame(me.predictor)
    n.pred<-length(pred[1,])
    pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred)
    me.pred<-matrix(data=me.pred[!is.na(me.pred)], ncol=n.pred)
    X<-cbind(1, pred)
    me1<-diag(me.response[!is.na(me.response)])
    x.ols<-cbind(1, pred)
    V1<-diag(rep(1, times=N))
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
    eY<-X%*%beta1;
    pred.mean<-X%*%beta1
    g.mean<-(t(rep(1, times=N))%*%solve(V1)%*%Y)/sum(solve(V1));
    sst<-t(Y-g.mean)%*% solve(V1)%*%(Y-g.mean)
    sse<-t(Y-pred.mean)%*%solve(V1)%*%(Y-pred.mean)
    sigma<-sse/(N-(n.pred+1))
    repeat{
      V<-diag(rep(sigma, times=N))+me1 + diag(as.numeric(me.pred%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),])))
      V.inverse<-solve(V)
      beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
      test<-matrix(nrow=(n.pred+1))
      for(f in 1:(n.pred+1))
      {
        if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
      }
      if(sum(test)==0) break
      beta1<-beta.i
    }
    beta1<-beta.i
    beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
    eY<-X%*%beta1;
    pred.mean<-X%*%beta1
    g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
    sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
    sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
    sigma<-sse/(N-(n.pred+1))
    r.squared<-(sst-sse)/sst

    sigma.sq.y<-(1/N)*(sum((Y-X%*%beta1)^2))
    log.like<-(-N/2)*log(2*pi)-(N/2)*log(sigma.sq.y)-sum((Y-X%*%beta1)^2/(2*sigma.sq.y))
    ml<-log.like
    reg.par<-matrix(data=0, nrow=(n.pred+1), ncol=2, dimnames=list(c("Intercept", if(n.pred==1) deparse(substitute(predictor)) else colnames(predictor)), c("Estimate", "Std. Error")))
    reg.par[,1]<-round(beta1,5);
    reg.par[,2]<-round(sqrt(diag(beta.i.var)),5)
    modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
    n.par=n.pred;
    modfit[1,1]=ml
    modfit[2,1]=-2*ml+2*(2+n.par)
    modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
    modfit[4,1]=-2*ml+log(N)*(2+n.par)
    modfit[5,1]=r.squared*100
    modfit[6,1]=sst
    modfit[7,1]=sse
    message("REGRESSION PARAMETERS");message("");
    print(reg.par);message("");
    message("--------------------------------------------------");
    message("MODEL FIT");message("");
    print(modfit); message("");
    message("==================================================");
  }

`oubm.sim` <-
  function(Xo, Yo, halflife, vx, vy, increment, b0, b1, topology, times, specs, n){
    alpha<-log(2)/halflife
    s.X<-sqrt(vx)
    s.y<-sqrt((vy*2*alpha))
    phylogeny<-cbind(topology, times, specs)
    N<-max(topology)  # N internal nodes
    int.node<-rep(NA, times=N)
    p<-matrix(data=0, nrow=length(specs), ncol=n)
    r<-matrix(data=0, nrow=length(specs), ncol=n)
    for(i in 1:n){
      z<-phylogeny.evolution(Xo, Yo, alpha, s.X, s.y, increment, b0, b1, topology, times)
      x.values<-z[,1]
      y.values<-z[,2]
      p[,i]<-as.numeric(c(int.node, x.values))
      r[,i]<-as.numeric(c(int.node, y.values))
    }
    sim.dat<-data.frame(cbind(phylogeny, p, r))
    x.head<-paste(rep("predictor_", times=n), 1:n, sep="")
    y.head<-paste(rep("response_", times=n), 1:n, sep="")
    names(sim.dat)<-c("ancestor", "time", "species", x.head, y.head)
    return(sim.dat)
  }

#' Title
#'
#' @param tree
#'
#' @return
#' @export
`ouch2slouch` <-
  function(tree)
  {

    if (!inherits(tree, "ouchtree"))
      stop(sQuote("tree"), " must be of class ", sQuote("ouchtree"))

    N<-length(tree@nodes)
    tmp<-as(tree, "data.frame")
    tmp$ancestors<-as.character(tmp$ancestors)
    tmp$ancestors<-as.numeric(tmp$ancestors)
    tmp$times<-as.character(tmp$times)
    tmp$times<-as.numeric(tmp$times)
    tmp$nodes<-as.character(tmp$nodes)
    tmp$nodes<-as.numeric(tmp$nodes)
    tmp$ancestors[1]<-0;
    slouch_node<-1:N
    ancestor<-rep(NA, times=N)

    for(i in 1:N)
    {
      if(tmp$labels[i]=="") tmp$labels[i]= NA;
    }

    rownames(tmp)<-1:nrow(tmp)
    names(tmp)=c("nodes", "ancestor", "time", "species")
    return(tmp)
  }

#' Parse tree function
#'
#' @param topology
#' @param times
#'
#' @return
#' @export
`parse.tree` <-
  function(topology, times){
    term<-terminal.twigs(topology);
    N <- length(term);
    dm<-distance.matrix(topology, times);
    bt<-branch.times(topology, times);
    pt<-list(N=N, term=term, dm=dm, bt=bt);
    return(pt);
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


`phylogeny.evolution` <-
  function(Xo, Yo, alpha, s.X, s.y, increment, b0, b1, topology, times){
    branch<-matrix(data=0, nrow=length(times), ncol=2);
    branch[1, 1]<-Xo;
    branch[1, 2]<-Yo;
    sx<-sqrt(increment)*s.X
    sy<-sqrt(increment)*s.y
    for(i in 2:length(times)){
      t<-times[i]-times[topology[i]];
      m<-round(t/increment)                   # note that some discretization error is introduced here, but it shouldnt matter as long as increment is small enough;
      X<-branch[topology[i], 1];
      Y<-branch[topology[i], 2];
      for (j in 1:m){
        Y<-Y+(-alpha*(Y-(b0+b1*(X)))*increment)+rnorm(1, mean=0, sd=sy);
        X<-X+rnorm(1, mean=0, sd=sx);
      }
      branch[i,1]<-X
      branch[i,2]<-Y
    }
    tip.specs=length(topology)-max(topology);
    branch1<-branch[tip.specs:length(topology),]
    return(branch1)
  }


pseudoinverse <-
  function (m){
    msvd <- svd(m)
    if (length(msvd$d) == 0) {
      return(array(0, dim(m)[2:1]))
    }
    else {
      return(msvd$v %*% (1/msvd$d * t(msvd$u)))
    }
  }

`regimes` <-
  function (topology, times, regime.specs, term) {
    N <- length(term);
    reg <- set.of.regimes(topology,regime.specs);
    R <- length(reg);
    beta <- vector(R*N, mode="list");
    for (i in 1:N) {
      for (k in 1:R) {
        p <- pedigree(topology, term[i]);
        n <- length(p);
        beta[[i + N*(k-1)]] <- as.integer(regime.specs[p[1:(n-1)]] == reg[k]);
      }
    }
    return(beta);
  }


`set.of.regimes` <-
  function (topology, regime.specs) {
    n <- length(regime.specs);
    id <- seq(1,n)[topology > 0];       # find all non-root nodes
    reg <- sort(unique(regime.specs[id]));
    return(reg);
  }

`sigma.X.estimate` <-
  function (predictor,me.predictor, topology, times) {
    pt <- parse.tree(topology,times);
    n <- pt$N;
    v <- pt$bt;
    w <- matrix(data=1,nrow=pt$N,ncol=1);
    me<-diag(me.predictor[!is.na(me.predictor)])
    dat <- predictor[!is.na(predictor)];
    beta<-solve(t(w)%*%solve(v)%*%w)%*%(t(w)%*%solve(v)%*%dat)
    e<-dat-beta
    sigma<-as.numeric((t(e)%*%solve(v)%*%e)/(n-1))
    repeat{
      beta<-solve(t(w)%*%solve(v + me/sigma)%*%w)%*%(t(w)%*%solve(v + me/sigma)%*%dat)
      e<-dat-beta
      theta <- beta
      sigma1<-(t(e)%*%solve(v +me/sigma)%*%e)/(n-1)
      if(abs(as.numeric(sigma1)-sigma) <= 0.0000001*sigma) break
      sigma <- as.numeric(sigma1)
    }
    return(list(as.numeric(theta), as.numeric(sigma)));
  }

#' Title
#'
#' @param topology
#' @param times
#' @param names
#' @param regimes
#' @param cex
#' @param lwd
#' @param reg.col
#'
#' @return
#' @export
slouchtree.plot <-function (topology, times, names = NULL, regimes = NULL, cex = NULL, lwd=NULL, reg.col=NULL) {
  if(is.null(cex)) cex<-1;
  if(is.null(lwd)) lwd<-1;
  rx <- range(times);
  rxd <- 0.1*diff(rx);

  if (is.null(regimes))
    regimes <- factor(rep(1,length(topology)));

  levs <- levels(as.factor(regimes));
  palette <- rainbow(length(levs));

  for (r in 1:length(levs)) {
    y <- tree.layout(topology);
    x <- times;
    f <- which(topology > 0 & regimes == levs[r]);
    pp <- topology[f];
    X <- array(data=c(x[f], x[pp], rep(NA,length(f))),dim=c(length(f),3));
    Y <- array(data=c(y[f], y[pp], rep(NA,length(f))),dim=c(length(f),3));
    oz <- array(data=1,dim=c(2,1));
    X <- kronecker(t(X),oz);
    Y <- kronecker(t(Y),oz);
    X <- X[2:length(X)];
    Y <- Y[1:(length(Y)-1)];
    if(!is.null(regimes))
    {if(is.null(reg.col))
      C <- rep(palette[r],length(X))
    }
    {if(!is.null(reg.col))
      C <- rep(reg.col[r],length(X))
    }
    if (r > 1) par(new=TRUE);
    par(yaxt='n')
    par(bty="n")
    par(font="2")

    plot(X,Y,type='l',col=C,lwd=lwd,xlab='time',ylab='',xlim = rx + c(-rxd,rxd),ylim=c(0,1));
    if (!is.null(names))
      text(X[seq(1,length(X),6)],Y[seq(1,length(Y),6)],names[f],pos=4, cex=cex);
  }
  par(yaxt="s") #reset graphic parameter to default
  par(bty="o")
  par(font="1")
}


#' Title
#'
#' @param topology
#'
#' @return
#' @export
`terminal.twigs` <-
  function(topology){
    n<-length(topology);
    return(seq(max(topology)+1, n));
  }

`tree.layout` <-
  function (topology) {
    root <- which(topology==0);
    return(arrange.tree(root,topology));
  }

#' Title
#'
#' @param topology
#' @param times
#'
#' @return
#' @export
`tsia` <-
  function (topology, times) {
    term <- terminal.twigs(topology);
    N <- length(term);
    t.ia <- matrix(data=0,nrow=N,ncol=N);
    t.ia[1,1]=0
    for (i in 2:N) {
      pedi <- pedigree(topology,term[i]);
      for (j in 1:(i-1)) {
        pedj <- pedigree(topology,term[j]);
        for (k in 1:length(pedi)) {
          if (any(pedj == pedi[k])) break;
        }
        t.ia[j,i] <- t.ia[i,j] <- (times[term[i]]-times[pedi[k]]);
      }
      t.ia[i,i] <- 0;
    }
    return(t.ia);
  }

#' Title
#'
#' @param topology
#' @param times
#'
#' @return
#' @export
`tsja` <-
  function (topology, times) {
    term <- terminal.twigs(topology);
    N <- length(term);
    t.ja <- matrix(data=0,nrow=N,ncol=N);
    t.ja[1,1]=0
    for (i in 2:N) {
      pedi <- pedigree(topology,term[i]);
      for (j in 1:(i-1)) {
        pedj <- pedigree(topology,term[j]);
        for (k in 1:length(pedi)) {
          if (any(pedj == pedi[k])) break;
        }
        t.ja[j,i] <- t.ja[i,j] <- (times[term[j]]-times[pedi[k]]);
      }
      t.ja[i,i] <- 0;
    }
    return(t.ja);
  }


make.slouch.data<-function(tree, data)
{
  if (inherits(tree, "phylo"))
  {
    tree<-ape2ouch(tree)
    tree<-ouch2slouch(tree)
  }
  if (inherits(tree, "ouchtree"))
  {
    tree<-ouch2slouch(tree)
  }

  N.col.tree<-length(tree[1,]);
  N.row.tree<-length(tree[,1]);
  N.col.data<-length(data[1,]);
  N.int<-length(tree$species[is.na(tree$species)]);
  int.tree<-tree[1:N.int,];
  tip.tree<-tree[(N.int+1):N.row.tree, ];
  tip.dat<-merge(tip.tree, data, by.x="species", by.y="species", sort=FALSE);
  x<-matrix(NA, ncol=(N.col.data-1), nrow=N.int);
  int<-data.frame(int.tree$species, int.tree$nodes, int.tree$ancestor, int.tree$time, x);
  names(int)<-names(tip.dat);
  slouch.dat<-rbind(int, tip.dat);

  return(slouch.dat)
}

is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))


rmNullObs <- function(x)
{
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}



####### WEIGHT MATRIX FOR FACTORS, BORROWED AND MODIFIED FROM OUCH #########
weight.matrix<-function(alpha, topology, times, N, regime.specs, fixed.cov, intercept)
{
  if (alpha == Inf) alpha<-10000000000000000000
  N <- N
  reg <- set.of.regimes(topology, regime.specs)
  R <- length(reg)
  T <- times[terminal.twigs(topology)]
  ep <- epochs(topology, times, terminal.twigs(topology))
  beta <- regimes(topology, times, regime.specs, terminal.twigs(topology))
  W <- matrix(data = 0, nrow = N, ncol = R + 1, dimnames = list(c(),
                                                                c("Ya", as.character(set.of.regimes(topology, regime.specs)))))
  W[, 1] <- exp(-alpha * T)
  for (i in 1:N) {
    delta <- diff(exp(alpha * (ep[[i]] - T[i])))
    for (k in 1:R) {
      W[i, k + 1] <- -sum(delta * beta[[i + N * (k - 1)]])
    }
  }
  if (is.null(intercept))
    W <- W
  else {
    if (intercept == "root") {
      root.reg <- as.character(regime.specs[times == 0])
      nonroot.reg <- as.character(reg[reg != root.reg])
      int <- as.matrix(W[, 1] + W[, root.reg])
      colnames(int) = root.reg
      if (max(int[, 1]) <= 0.01)
        W <- W
      else {
        W2 <- cbind(int, W[, nonroot.reg])
        W <- W2
        if(max(abs(W[,1]-W[,2]))<= 0.01) W<-W[,(2:length(W[1,]))] #sjekk denne for treghet

        #W <- W2[, 2:length(W2[1,])]
      }
    }
    else W[, 1] <- intercept
  }
  if (max(W[, 1]) <= 0.01)
    W <- W[, -1];

  if(!is.null(fixed.cov))
  {
    fixed.pred<-data.frame(fixed.cov);
    n.fixed.pred<-length(fixed.pred[1,]);
    fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred);
    W<-cbind(W, fixed.pred)
  }
  return(W)
}

np.regression<-function(response, me.response, predictor, me.predictor, convergence=NULL){
  if(is.null(convergence)) convergence=0.000001
  Y <- response[!is.na(response)];
  N <- length(Y);
  pred<-data.frame(predictor)
  me.pred<-data.frame(me.predictor)
  n.pred<-length(pred[1,])
  pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred)
  me.pred<-matrix(data=me.pred[!is.na(me.pred)], ncol=n.pred)
  X<-cbind(1, pred)
  me1<-diag(me.response[!is.na(me.response)])
  x.ols<-cbind(1, pred)
  V1<-diag(rep(1, times=N))
  beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  eY<-X%*%beta1;
  pred.mean<-X%*%beta1
  g.mean<-(t(rep(1, times=N))%*%solve(V1)%*%Y)/sum(solve(V1));
  sst<-t(Y-g.mean)%*% solve(V1)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V1)%*%(Y-pred.mean)
  sigma<-sse/(N-(n.pred+1))
  repeat{
    V<-diag(rep(sigma, times=N))+me1 + diag(as.numeric(me.pred%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),])))
    V.inverse<-solve(V)
    beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
    test<-matrix(nrow=(n.pred+1))
    for(f in 1:(n.pred+1))
    {
      if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
    }
    if(sum(test)==0) break
    beta1<-beta.i
  }
  beta1<-beta.i
  beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
  eY<-X%*%beta1;
  pred.mean<-X%*%beta1
  g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
  sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
  sigma<-sse/(N-(n.pred+1))
  r.squared<-(sst-sse)/sst

  sigma.sq.y<-(1/N)*(sum((Y-X%*%beta1)^2))
  log.like<-(-N/2)*log(2*pi)-(N/2)*log(sigma.sq.y)-sum((Y-X%*%beta1)^2/(2*sigma.sq.y))
  ml<-log.like
  reg.par<-matrix(data=0, nrow=(n.pred+1), ncol=2, dimnames=list(c("Intercept", if(n.pred==1) deparse(substitute(predictor)) else colnames(predictor)), c("Estimate", "Std. Error")))
  reg.par[,1]<-round(beta1,5);
  reg.par[,2]<-round(sqrt(diag(beta.i.var)),5)
  modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
  n.par=n.pred;
  modfit[1,1]=ml
  modfit[2,1]=-2*ml+2*(2+n.par)
  modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
  modfit[4,1]=-2*ml+log(N)*(2+n.par)
  modfit[5,1]=r.squared*100
  modfit[6,1]=sst
  modfit[7,1]=sse
  message("REGRESSION PARAMETERS");message("");
  print(reg.par);message("");
  message("--------------------------------------------------");
  message("MODEL FIT");message("");
  print(modfit); message("");
  message("==================================================");
}




