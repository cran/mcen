#' Provides initial estimates for the mcen functionF
#'@param x the n x p design matrix
#'@param y the n x y matrix of responses
#'@param family type of likelihood used two options "mgaussian" or "mbinomial"
#'@param delta sparsity tuning parameter
#'@param gamma_y tuning parameter for clustering responses
#'@param intercept whether an intercept should be included in the model
#'@return matrix of coefficients
#'@importFrom glmnet glmnet
#'@importFrom stats coefficients
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu>        
mcen.init <- function(x,y,family="mgaussian",delta=NULL,gamma_y=1,intercept=FALSE){
    if(length(delta) > 1 | length(gamma_y) > 1){
        stop('initialization only for single delta, gamma_y')
    }
    r <- dim(y)[2]
    n <- dim(y)[1]
    p <- dim(x)[2]
    if(family=="mgaussian"){
        init_family = "gaussian"
        intercept_b=intercept
    }
    if(family=="mbinomial"){
      init_family="binomial"
      intercept_b=TRUE
      x=x[,-1]
      ## probably need to force the intercept here to?
    }
    init_model <- apply(y,2,glmnet,family=init_family,x=x,alpha=delta/(delta+gamma_y),lambda=(delta+gamma_y)/(2*n),intercept=intercept_b)
    if(intercept_b){
        do.call(cbind,lapply(init_model,coefficients))
    } else{
        do.call(cbind,lapply(init_model,coefficients))[-1,]
    }
}

#' Fits an MCEN model 
#'@param x Matrix of predictors.
#'@param y Matrix of responses. 
#'@param family Type of likelihood used two options "mgaussian" or "mbinomial".
#'@param ky Clusters for response.
#'@param delta L1 penalty.
#'@param gamma_y Penalty for with y clusters difference in predicted values.
#'@param ndelta Number of delta parameters.
#'@param delta.min.ratio Ratio between smallest and largest delta.
#'@param eps Convergence criteria.
#'@param scale_x Whether x matrix should be scaled, default is True. 
#'@param scale_y Whether y matrix should be scaled, default is True. 
#'@param clusterMethod K-means function used kmeans or kmeanspp.
#'@param clusterStartNum Number of random starting points for clustering.
#'@param clusterIterations Number of iterations for cluster convergence.
#'@param cluster_y An a priori definition of clusters. If clusters are provided they will remain fixed and are not estimated. Objective function is then convex. 
#'@param max_iter Maximum number of iterations for coefficient estimates.
#'@param init_beta Clustering step requires an initial estimate, default is to use elastic net solution.
#'@param n.cores Number of cores used for calculation default is 1.
#'@return returns a MCEN object
#' \item{beta}{List of the coefficient estimates.}
#' \item{delta}{Value of delta.}
#' \item{gamma_y}{Value of gamma_y.}
#' \item{ky}{Value of ky.} 
#' \item{y_clusters}{List of the clusters of y.} 
#'@examples
#' x <- matrix(rnorm(400),ncol=4)
#' beta <- beta <- matrix(c(1,1,0,0,0,0,-1,-1,0,0,-1,-1,1,1,0,0),ncol=4)
#' y <- x%*%beta + rnorm(400) 
#' mcen_fit <- mcen(x,y,ky=2,delta=1)
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
#'@references Price, B.S. and Sherwood, B. (2018). A Cluster Elastic Net for Multivariate Regression. arXiv preprint arXiv:1707.03530. \url{http://arxiv-export-lb.library.cornell.edu/abs/1707.03530}.
#'@export
#'@importFrom stats cov

mcen <- function(x,y,family="mgaussian",ky=NULL,delta=NULL,gamma_y=1,ndelta=25,delta.min.ratio = NULL,
						eps=.00001,scale_x=TRUE, scale_y=TRUE, clusterMethod="kmeans",clusterStartNum=30,clusterIterations=10,cluster_y=NULL,
						max_iter=10,init_beta=NULL, n.cores=1){  
  if(is.null(cluster_y)==FALSE){
	cluster_y_num <- length(unique(cluster_y))
	if(is.null(ky)==FALSE){
		if(cluster_y_num !=ky){
			ky <- cluster_y_num
			warning("ky does not match number of clusters in cluster_y, cluster_y then overrides ky")
		}
	}
  }
  if(scale_x){
    x <- scale(x)
    sigma_x <- attributes(x)$'scaled:scale'
    mean_x <- attributes(x)$'scaled:center'
  }
  if(scale_y && family!="mbinomial"){
    y <- scale(y)
    mean_y <- attributes(y)$'scaled:center'
    sigma_y  <- attributes(y)$'scaled:scale'
  }
  p=dim(x)[2]
  r=dim(y)[2]
  n=dim(x)[1]
   if(is.null(delta.min.ratio)==TRUE){
    delta.min.ratio=ifelse(n < p, 0.01, 1e-04)
  }
  if(family=="mgaussian"){
  Cxy=cov(x,y)*(n-1)
  Cxx=cov(x)*(n-1)
  }
  if(is.null(delta)){
  if(family=="mgaussian"){
  delta_max = 2*max(abs(Cxy))+1
  }
  if(family=="mbinomial"){
      ## Need to add bens rule here
      phat<-apply(y,2,mean)
      zw<-log(phat/(1-phat))*phat*(1-phat)+y-phat
      delta_max=2*max(abs(apply(zw,2,function(m){t(x)%*%zw})))+1
    }
      # plus one is just to make sure we get a full sparse matrix returned
	delta_min = delta.min.ratio*delta_max
	delta = exp(seq(log(delta_max),log(delta_min),length.out=ndelta))
  }
  
  ndelta <- length(delta)
  
  delta_iter <- 1
  beta_list <- list()
  cluster_list <- list()
  init_beta <- NULL
  if(family=="mbinomial"){
    x=cbind(rep(1,n),x)
  }	
  for(i in 1:ndelta){
	if(i == 1 || sum(init_beta)==0){
  	#first step get initial estimates, however we don't want to get stuck at initial estimates of all zeros for coefficients
	  if(length(delta)>1){
	  init_beta <- mcen.init(x,y,family,delta[i],gamma_y)
	  }else{
	    init_beta<-mcen.init(x,y,family,delta,gamma_y)
	  }
  	}
  	if(family=="mgaussian"){
    beta_cluster <- mcen_workhorse(init_beta,delta[i],Cxx,Cxy,family,ky=ky,gamma_y=gamma_y,eps=eps,
  										clusterMethod=clusterMethod,clusterIterations=clusterIterations,
  										clusterStartNum=clusterStartNum,cluster_y=cluster_y, max_iter=max_iter, x=x)
  	}
    if(family=="mbinomial"){
      
      if(length(delta)>1){
      beta_cluster<-mcen_bin_workhorse(init_beta,delta=delta[i],x=x,y=y,family="mbinomial",ky=ky,gamma_y=gamma_y,eps=eps,
                                      clusterMethod=clusterMethod,clusterIterations=clusterIterations,
                                       clusterStartNum=clusterStartNum,cluster_y=cluster_y, max_iter=max_iter)
      }else{
        beta_cluster<-mcen_bin_workhorse(init_beta,delta=delta,x=x,y=y,family="mbinomial",ky=ky,gamma_y=gamma_y,eps=eps,
                                         clusterMethod=clusterMethod,clusterIterations=clusterIterations,
                                         clusterStartNum=clusterStartNum,cluster_y=cluster_y, max_iter=max_iter)
      }
    }

    init_beta <- beta_cluster$beta
    if(family=="mgaussian"){
  	beta_list[[i]] <- beta_adjust(init_beta,sigma_x,sigma_y,mean_x,mean_y)
    }
    if(family=="mbinomial"){
      beta_list[[i]]<-beta_adjust_bin(init_beta,sigma_x)
    }
  	cluster_list[[i]] <- beta_cluster$clusters
  }
  return_val <- list(beta=beta_list,delta=delta,gamma_y=gamma_y,ky=ky,y_clusters=cluster_list)
  if(family=="mgaussian"){
  class(return_val) <- c("mgauss_mcen", "mcen")
  }
  if(family=="mbinomial"){
    class(return_val)<-c("mbinom_mcen","mcen")
  }
  return_val
}

#' Estimates the clusters and provides the coefficients for an mcen object 
#'@param beta The initial value of the coefficients
#'@param delta The sparsity (L1) tuning parameter
#'@param xx Matrix of transpose of x times x. 
#'@param xy Matrix of transpose of x times y. 
#'@param family Type of likelihood used two options "mgaussian" or "mbinomial" 
#'@param ky Number of clusters for the response 
#'@param gamma_y Penalty for the y clusters difference in predicted values
#'@param eps Convergence criteria
#'@param clusterMethod Which clustering method was used, currently support kmeans or kmeanspp
#'@param clusterIterations Number of iterations for cluster convergence
#'@param clusterStartNum Number of random starting points for clustering
#'@param cluster_y An a priori definition of clusters. If clusters are provided they will remain fixed and are not estimated. Objective function is then convex. 
#'@param max_iter The maximum number of iterations for estimating the coefficients
#'@param x The design matrix    
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
#'@importFrom methods as
#'@importFrom Matrix sparseMatrix
#'@useDynLib mcen CDU
mcen_workhorse<-function(beta,delta=NULL,xx,xy,family="mgaussian",ky=NULL,gamma_y=.5,eps=.00001,
                            clusterMethod="kmeans",clusterIterations=100,clusterStartNum=30,cluster_y=NULL,max_iter=10, x=x){
  #really need to think about how to do this function if x and y are not scaled and there is an intercept. This maybe important for the glm case
  beta_dim <- dim(beta)
  p <- beta_dim[1]
  r <- beta_dim[2]
  if(sum(beta==0)==p*r){
    #if initial estimate is totally sparse then mcen update is not needed. Philosophically, I'm ok with this because mcen idea is 
    # to group coefficient estimates, if all coefficients are zero there is nothing to group and nothing for mcen to do. 
    return_val <- list(beta = as(beta,"sparseMatrix"), clusters=rep(0,r))
	return( return_val)
  } else{
    fitted_vals <- matrix(x%*%beta,ncol=r)
    iter <- 0
    fitting_done <- FALSE
    while(!fitting_done){
      #print(paste("fitting done is ", fitting_done, " and iter is ", iter))
	  num_potential_centers <- ncol( unique(fitted_vals,MARGIN=2))
	  if(num_potential_centers < ky){
		fitting_done <- TRUE
		warning(paste0("Algorithm stopped early for gamma_y ", gamma_y, ", delta ", delta, " and ky ", ky, " because number of unqiue fitted values is less than ", ky, " Common for large values of deltas that cause spares solution and little variability in fitted values.\n")) 
		clusters_used <- rep(0,r)
		y_clusters <- NULL
	  } else if(is.null(cluster_y)==FALSE){
		clusters_used <- cluster_y
	  } else{
		y_clusters <- cluster(fitted_vals,ky,clusterMethod,clusterIterations,clusterStartNum)
		clusters_used <- y_clusters$cluster
	  }
      if(iter > 0 && fitting_done==FALSE){
		if(is.null(cluster_y)){
			fitting_done=SetEq(y_clusters2$cluster,y_clusters$cluster)
			#fitting_done=SetEq(y_clusters2$centers,y_clusters$centers)
		} else{
			fitting_done <- TRUE
		}
	  }
	  if(fitting_done == FALSE){
		  beta <-matrix(.C("CDU",beta=as.double(as.vector(beta)),xx=as.double(as.vector(xx)), xy=as.double(as.vector(xy)), y_clusters=as.double(clusters_used),#y_clusters$cluster), 
							delta=as.double(delta), gamma_y=as.double(gamma_y), eps=as.double(eps),miter=as.integer(max_iter), r=as.integer(r), p=as.integer(p), beta0=as.double(as.vector(beta+1)), set=as.double(rep(0,r)), First=as.double(0),
						    Second=as.double(0),mine=as.double(rep(0,p*r)))$beta, nrow=p,ncol=r,byrow=FALSE)
	  } 
	  if(iter == max_iter){
			fitting_done <- TRUE
	  }
	  iter <- iter + 1
	  if(is.null(cluster_y)==TRUE){
		y_clusters2 <- y_clusters
	  }
	  fitted_vals=matrix(x%*%beta,ncol=r)
	}
	   if(is.null(cluster_y)==FALSE){
			fitting_done <- TRUE #only one iteration of fitting done if clusters are prespecified
	   }
    }
	return_val <- list(beta = as(beta,"sparseMatrix"), clusters=clusters_used)
    return( return_val)
}							

#' Adjusts the value of the coefficients to account for the scaling of x and y.
#'@param beta The estiamte of beta with scaled data.
#'@param sigma_x Sample tandard deviations of the original predictors.
#'@param sigma_y Sample standard deviations of the orignal responses.
#'@param mean_x Sample means of the original predictors .
#'@param mean_y Sample means of the original responses.
#'@return Returns the adjusted coefficients 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
beta_adjust <- function(beta,sigma_x,sigma_y,mean_x,mean_y){
	beta <- beta*sigma_x
	for(i in 1:dim(beta)[2]){
      beta[,i] <- beta[,i]*sigma_y[i]
    }
	intercept <- mean_y - mean_x %*% beta
	beta <- rbind(intercept,beta)
	beta
}

#' Adjusts the value of the binomial coefficients to account for the scaling of x.
#'@param beta The estiamte of beta with scaled data.
#'@param sigma_x Sample tandard deviations of the original predictors.
#'@return Returns the adjusted coefficients 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
beta_adjust_bin<-function(beta,sigma_x){
  beta_out<-beta
  beta_out[-1,]-beta[-1,]*sigma_x
  return(beta_out)
}

#' SetEq test set equivalence of two clustering sets
#'@param set1 is the cluster assignments of the previous iteration
#'@param set2 is the cluster assignments of the current clusters
#'@return Returns a logical saying if the two clusterings are equal
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
SetEq <-function(set1,set2){
  s1<-vector("list",max(set1))
  for(i in 1:max(set1)){
    s1[[i]]=which(set1==i)
  }
  Sets=sapply(s1,function(z){
    s2=vector("list",max(set2))
    for(i in 1:max(set2)){
      s2[[i]]=which(set2==i)
    }
    sum(sapply(s2,function(w){identical(z,w)}))
  })
  sum(Sets)==max(set1)
}



#' Wrapper function for different clustering methods
#'@param x data to be clustered. Clustering will be done on the columns.
#'@param cNum number of cluster centers
#'@param clusterMethod "kmean" for kmeans function, "kmeanspp" for kcca implementation of kmeans++
#'@param clusterIterations number of maximum iterations for clustering
#'@param clusterStartNum random number of starting points used
#'@return Returns cluster assignments
#'@importFrom flexclust kcca
#'@importFrom stats kmeans
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
cluster <- function(x, cNum, clusterMethod="kmeans", clusterIterations=100, clusterStartNum=30){
#get cluster results using method on data x using cluster tuning paramater c_tune
    if(clusterMethod == "kmeans"){
      kmx <- kmeans(t(x), cNum, iter.max=clusterIterations, nstart=clusterStartNum)
      return_val <- list(cluster=kmx$cluster, centers=t(kmx$centers), av_dist=(kmx$withinss/kmx$size)/dim(x)[1], size=kmx$size)
    }
    else if(clusterMethod == "kmeanspp"){                                                         
       kmx <- kcca(t(x),cNum,control=list(initcent="kmeanspp",iter.max=clusterIterations))
       centers <- t(attributes(kmx)$centers)
       #center_trans <- centers[,attributes(kmx)$cluster]
       #av_dist <- tapply(colMeans((x-center_trans)^2),attributes(kmx)$cluster,mean)
       return_val <- attributes(kmx)$cluster
    }
    return_val
}

#' randomly assign n samples to k groups
#' @param n number of samples
#' @param k number of groups
#' @return Returns assignments of n into k groups
#' @author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
randomly_assign <- function(n,k){
#randomly assign n samples into k groups
   small_set <- floor(n/k)
   group_assign <- NULL
  if(n %% k == 0){
     group_assign <-  rep(seq(1,k),n/k)
   } else{
     remainder <- n %% k
     for(i in 1:remainder){
        group_assign <- c(group_assign, rep(i,small_set+1))
     }
     group_assign <- c(group_assign, rep(seq((i+1),k),small_set))
   }
   sample(group_assign)
}

#'matrix multiply
#'@param beta Matrix of coefficients.
#'@param x Design matrix. 
#'@return Returns x times beta
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
matrix_multiply <- function(beta,x){
	x%*%beta
}

#' predictions from a mcen model
#'@method predict mcen
#'@param object The mcen object.
#'@param newx A matrix of new observations.
#'@param ... Additional variables to be sent to predict. 
#'@return Returns predictions for each beta of an mcen object
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu>
#'@examples
#' x <- matrix(rnorm(400),ncol=4)
#' beta <- beta <- matrix(c(1,1,0,0,0,0,-1,-1,0,0,-1,-1,1,1,0,0),ncol=4)
#' y <- x%*%beta + rnorm(400) 
#' mcen_fit <- mcen(x,y,ky=2,delta=1)
#' new_x <- matrix(rnorm(12),ncol=4)
#' mcen_preds <- predict(mcen_fit, new_x)
#'@export
predict.mcen <- function(object,newx,...){
	newx <- cbind(1,newx)
	preds <- lapply(object$beta,matrix_multiply,newx)
	preds
}

#' Calculates the out of sample likelihood for an mcen object
#' @param obj The mcen object.
#' @param test_x The matrix of test predictors.
#' @param test_y The matrix of test responses. 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
pred_eval <- function(obj,test_x,test_y) 
{
    UseMethod("pred_eval")
}

#' Calculates sum of squared error between two vectors or matrices
#' @param pred the predictions
#' @param test_y the testing response values
#' @return returns the sum of the squared differences between pred and test_y
#' @author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
squared_error <- function(pred,test_y){
	sum((pred-test_y)^2)
}

#' Calculates the prediction error for a mgauss_mcen object. 
#'@param obj The mgauss_mcen object.
#'@param test_x The matrix of test predictors. 
#'@param test_y The matrix of test responses. 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
#'@importFrom stats predict
pred_eval.mgauss_mcen <- function(obj,test_x,test_y){
	predicted_vals <- predict(obj,test_x)
	pred_error <- sapply(predicted_vals,squared_error,test_y)
	pred_error 
}

#' Calculates out of sample error on the binomial likelihood
#'@param pred The predicted values. 
#'@param test_y The test response values.  
#'@importFrom faraway ilogit
#'@author Brad Price <brad.price@@mail.wvu.edu> 
vl_binom<-function(pred,test_y){
    pred=as.matrix(pred)
    out=which(ilogit(pred)=="NaN")
    pred[out]=1-10^-10
    pred[which(ilogit(pred)>1-10^-10)]=1-10^-10
    pred[which(ilogit(pred)<10^-10)]=10^-10
    foo=sum(log(ilogit(pred)^test_y)+log((1-ilogit(pred))^(1-test_y)))
    return(-foo)
}

#' Evaluates prediction error for multiple binomial responses.
#'@param obj The mbinom_mcen object. 
#'@param test_x A matrix of the test predictors. 
#'@param test_y A matrix of the test responses. 
#' @author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
#'@importFrom stats predict
pred_eval.mbinom_mcen<-function(obj,test_x,test_y){
  predicted_vals<-predict(obj,test_x)
  vl<-sapply(predicted_vals,vl_binom,test_y)
  return(vl)
}

#' Cross validation for mcen function
#' @param x Matrix set of predictors.
#' @param y Matrix set of responses.
#' @param family The exponential family the response corresponds to.
#' @param ky A vector with the number of possible clusters for y.
#' @param gamma_y Set of tuning parameter for clustering penalty in response categories.
#' @param nfolds Number of folds used in the cross-validation. 
#' @param folds A vector of length n, where this identifies what fold of the kfold cross validation each observation belongs to.
#' @param cluster_y a priori definition of clusters. If clusters are provided they will remain fixed and are not estimated. Objective function is then convex. 
#' @param n.cores Number of cores used for parallel processing. 
#'@param ... The variables passed to mcen
#' @return Returns a cv.mcen object. 
#' \item{models}{A list of mcen objects.}
#' \item{cv}{Cross validation results.} 
#' \item{ky}{The same value as the input ky.} 
#' \item{gamma_y}{The same value as the input gamma_y.} 
#'@examples
#'\donttest{
#' x <- matrix(rnorm(400),ncol=4)
#' beta <- beta <- matrix(c(1,1,0,0,0,0,-1,-1,0,0,-1,-1,1,1,0,0),ncol=4)
#' y <- x%*%beta + rnorm(400) 
#' cv_fit <- cv.mcen(x,y,ky=2)
#'}
#'@references Price, B.S. and Sherwood, B. (2018). A Cluster Elastic Net for Multivariate Regression. arXiv preprint arXiv:1707.03530. \url{http://arxiv-export-lb.library.cornell.edu/abs/1707.03530}.
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
#'@export
#'@importFrom parallel mclapply
cv.mcen <-function(x,y,family="mgaussian",ky=seq(2,4),gamma_y=seq(.1,5.1,.5),nfolds=10,folds=NULL,cluster_y=NULL,delta=NULL, n.cores=1,...){

  
  if(is.null(cluster_y)==FALSE){
	ky = length(unique(cluster_y))
	warning("ky does not match number of clusters in cluster_y, cluster_y then overrides ky")
  }
  clus=length(ky)
  gamma_yl=length(gamma_y)
  
  ResOut <- NULL#vector("list",clus)
  
  if(is.null(folds)){
	folds <- randomly_assign(dim(x)[1],nfolds)
  }
  
  Eval<-function(del,gammay,kyn,fam,...){
	 #print(paste("working on eval for delta", del, "and gamma", gammay))
     pse <- 0
	 for(s in unique(folds)){
		 #print(paste("working on fold", s))
		 Set=which(folds==s)
		 #bet=mcen(x[-Set,],y[-Set,],delta=del,gamma_y=gammay,ky=kyn,...)
		 if(fam=="mbinomial"){
		   train_model <- mcen(x[-Set,],y[-Set,],family="mbinomial",delta=del,gamma_y=gammay,ky=kyn,cluster_y=cluster_y,...)
		 }else{
		   train_model <- mcen(x[-Set,],y[-Set,],delta=del,gamma_y=gammay,ky=kyn,cluster_y=cluster_y,...) 
		 }
		 pse <- pse + pred_eval(train_model,x[Set,],y[Set,])
	 }
	 pse
  }
  
  mc_Eval <- function(mcen_obj,family,...){
  	delta_vals <- mcen_obj$delta
  	gamma_y <- mcen_obj$gamma_y
  	ky <- mcen_obj$ky
  	#print(paste("Evaluating model of ky working on ", ky, " gammay is ", gamma_y))
  	pe <- Eval(delta_vals,gamma_y,ky,fam=family,...)
    cbind(ky,gamma_y,delta_vals,pe)
  }
  
  mc_mcen <- function(tune_vals,family,...){
  	#print(paste("Fitting full model ky",tune_vals[1],"gammay",tune_vals[2]))
  	mcen(x,y,ky=tune_vals[1],gamma_y=tune_vals[2],cluster_y=cluster_y,family=family,...)
  }
  
  mcen_list <- list()
	mcen_num <- 1
	if(n.cores == 1){ 
		for(k in 1:clus){
			k_gamma_results <- NULL
			for(g in 1:gamma_yl){
				#print(paste0("working on k ", k, " and g ", g))
				mcen_full <- mcen(x,y,family,ky=ky[k], gamma_y=gamma_y[g],cluster_y=cluster_y,delta=delta,...)
				mcen_list[[mcen_num]] <- mcen_full
				train_delta <- mcen_full$delta
				eval_results <- Eval(train_delta,gamma_y[g],ky[k],fam=family,...)
				ResOut <- rbind(ResOut, cbind(ky[k],gamma_y[g],train_delta,eval_results,mcen_num))
				mcen_num <- mcen_num + 1
			}
		} 
	}  else{
		ky_vals <- sort(rep(ky,gamma_yl))
		gamma_y_vals <- rep(gamma_y,clus)
		tune_mat <- rbind(ky_vals,gamma_y_vals)
		tune_list <- split(tune_mat, rep(1:ncol(tune_mat), each = nrow(tune_mat)))
		
		mcen_list <- mclapply(tune_list,mc_mcen,mc.cores=n.cores,family=family,...)
		ResOut <- mclapply(mcen_list,mc_Eval,mc.cores=n.cores,family=family,...)
		ResOut <- do.call(rbind,ResOut)
		model_num <- rep(seq(1,length(mcen_list)),each=length(mcen_list[[1]]$delta))
		ResOut <- cbind(ResOut,model_num)
    
	 }
 
 
  #returns list of all models and cross-validation erorr
  #later functions will handle coefficients and predictions from this object
  Out = list(models=mcen_list,cv=ResOut,ky=ky,gamma_y=gamma_y)
  class(Out) <- c(paste0(family,"_cv.mcen"),"cv.mcen")
  return(Out)

}


#' Prints nice output for an mcen object. 
#'@param x The mcen object. 
#'@param ... Additional parameters. 
#'@return Prints out some basic information about the mcen object. 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
#'@export
print.mcen <- function(x,...){
  cat("\n mcen object \n")
  cat("class:")
  print(class(x))
  cat("\n gamma_y value")
  print(x$gamma_y)
  cat("\n Number of groups for response:")
  print(x$ky)
}

#' Prints nice output for a cv.mcen object. 
#'@param x The cv.mcen object. 
#'@param ... Additional parameters.
#'@return Prints out information about where the cv.mcen object was minimized. 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu> 
#'@export
print.cv.mcen <- function(x,...){
   cat("\nCross Validation Error minimized at:\n")
   min_spot <- which.min(x$cv[,4])
   print(x$cv[min_spot,])
}

#' Makes predictions from the model with the smallest cross-validation error. 
#'@method predict cv.mcen
#'@param object The cv.mcen object. 
#'@param newx The X matrix of predictors. 
#'@param ... Additional parameters to be sent to predict. 
#'@return Returns the predicted values from the model with the smallest cross-validation error. 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu>
#'@examples
#'\donttest{
#' x <- matrix(rnorm(400),ncol=4)
#' beta <- beta <- matrix(c(1,1,0,0,0,0,-1,-1,0,0,-1,-1,1,1,0,0),ncol=4)
#' y <- x%*%beta + rnorm(400) 
#' mcen_fit <- cv.mcen(x,y,ky=2,gamma_y=3)
#' new_x <- matrix(rnorm(12),ncol=4)
#' mcen_preds <- predict(mcen_fit, new_x)
#'}
#'@export
#'@importFrom stats coef
predict.cv.mcen <- function(object,newx,...){
  newx <- cbind(1,newx)
  preds <- newx %*% coef(object)
  if("mbinomial_cv.mcen" %in% class(object)){
    ilogit(preds)
  } else{
    preds
  }
}

#' Returns the cluster values from a cv.mcen object. 
#'@param obj The cv.mcen object. 
#'@return Returns the clusters from the model with the smallest cross-validation error. 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu>
#'@examples
#'\donttest{
#' x <- matrix(rnorm(400),ncol=4)
#' beta <- beta <- matrix(c(1,1,0,0,0,0,-1,-1,0,0,-1,-1,1,1,0,0),ncol=4)
#' y <- x%*%beta + rnorm(400) 
#' mcen_fit <- cv.mcen(x,y,ky=2,gamma_y=3)
#' mcen_cluster <- cluster.vals(mcen_fit)
#' }
#'@export
cluster.vals <- function(obj){
   min_spot <- which.min(obj$cv[,4])
   min_delta <- obj$cv[min_spot,3]
   min_model <- obj$cv[min_spot,5]
   best_model <- obj$models[[min_model]]
   delta_spot <- which(best_model$delta==min_delta)
   best_model$y_clusters[[delta_spot]]
}

#' Returns the coefficients from the cv.mcen object with the smallest cross-validation error.
#'@method coef cv.mcen
#'@param object The cv.mcen object.
#'@param ... Additional values to be passed. 
#'@return The matrix of coefficients for the best MCEN model as determined by cross-validation. 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu>
#'@examples
#'\donttest{
#' x <- matrix(rnorm(400),ncol=4)
#' beta <- beta <- matrix(c(1,1,0,0,0,0,-1,-1,0,0,-1,-1,1,1,0,0),ncol=4)
#' y <- x%*%beta + rnorm(400) 
#' mcen_fit <- cv.mcen(x,y,ky=2,gamma_y=3)
#' best_coef <- coefficients(mcen_fit)
#' }
#'@export
#'@importFrom stats coef
coef.cv.mcen <- function(object,...){
   min_spot <- which.min(object$cv[,4])
   min_delta <- object$cv[min_spot,3]
   min_model <- object$cv[min_spot,5]
   best_model <- object$models[[min_model]]
   best_beta <- coef(best_model,min_delta) 
   best_beta
}

#' Returns the coefficients from an mcen object. 
#'@method coef mcen 
#'@param object The mcen object. 
#'@param delta The L1 tuning parameter 
#'@param ... Additional values to pass on. 
#'@return The matrix of coefficients. 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu>
#'@examples 
#' x <- matrix(rnorm(400),ncol=4)
#' beta <- beta <- matrix(c(1,1,0,0,0,0,-1,-1,0,0,-1,-1,1,1,0,0),ncol=4)
#' y <- x%*%beta + rnorm(400) 
#' mcen_fit <- mcen(x,y,ky=2,gamma_y=3,delta=c(1,2))
#' best_coef <- coefficients(mcen_fit,delta=1)
#'@export
coef.mcen <- function(object,delta=NULL,...){
	if(is.null(delta)){
		object$beta
	} else{
		delta_pos <- which(object$delta == delta)
		object$beta[[delta_pos]]
	}
}

#' Gets the index position for the model with the smallest cross-validation error. 
#'@param model The cv.mcen object.
#'@return Returns the index for the model with the smallest cross-validation error.
#'@author Ben Sherwood <ben.sherwood@@ku.edu>, Brad Price <brad.price@@mail.wvu.edu>
#'@examples
#'\donttest{
#' x <- matrix(rnorm(400),ncol=4)
#' beta <- beta <- matrix(c(1,1,0,0,0,0,-1,-1,0,0,-1,-1,1,1,0,0),ncol=4)
#' y <- x%*%beta + rnorm(400) 
#' mcen_fit <- cv.mcen(x,y,ky=2,gamma_y=3)
#' get_best_cvm(mcen_fit)
#' }
#'@export
get_best_cvm <- function(model){
	min(model$cvm)
}




#'Creates the probabilities and working response for the glmnet update for a given response with a binomial family
#' @param X the matrix of predictors. 
#' @param Beta current iteration of the regression coefficients
#' @param Y is the matrix of responses
#' @param r the response of interest  
#' result is a list of things needed for the working response in glmnet
#'@author Brad Price <brad.price@@mail.wvu.edu>
CalcHorseEBin<-function(X,Beta,Y,r){
  Probs=ilogit(X%*%Beta[,r])
  set=which(Probs<=10^-5 | Probs>=(1-10^-5))
  W=Probs*(1-Probs)
  W[set]=10^-5
  Z=X%*%Beta[,r]+(Y[,r]-Probs)/W
  V=sqrt(W)*Z
  Tm=apply(X,2,function(x){sqrt(W)*x})
  R=t(X)%*%X
  return(list(P=Probs,W=W,Z=Z,V=V,Tm=Tm,R=R))
}


#' Creates the the working response for all responses for glmnet binomial family
#' @param X the matrix of predictors.
#' @param Beta current iteration of the regression coefficients
#' @param Y is the matrix of responses
#' result is the list of vectors needed for the working responses in glmnet
#'@author Brad Price <brad.price@@mail.wvu.edu>
CalcHorseBin<-function(Y,X,Beta){
  TMm=NULL
  TmV=NULL
  Cm=NULL
  Vm=NULL
  R=NULL
  for(l in 1:dim(Y)[2]){
    out<-CalcHorseEBin(X,Beta,Y,l)
    TMm=cbind(TMm,out$Tm)
    Cm=cbind(Cm,t(out$Tm)%*%out$Tm)
    Vm=cbind(Vm,out$V)
    TmV=cbind(TmV,t(out$Tm)%*%out$V)
    ##R=cbind(R,out$R)
  }
  R=out$R
  return(list(Tm=TMm,C=Cm,V=Vm,R=R,TmV=TmV))
}


#' The workhorse function for the binomial updates in mcen.  It uses IRWLS glmnet updates to solve the regression problem.
#' @param Y the matrix of responses
#' @param X the matrix of predictors with the intercept included
#' @param delta the tuning parameter for the lasso penalty
#' @param gamma_y the tuning parameter for the ridge fusion penalty
#' @param y_clusters the cluster assignments from the provided clustering algorithm 
#' @param set_length the size of each cluster corresponding to a given response.  r dimensions with each element containing the cluster size of that responses cluster.
#' @param eps the tolerance for conversion normally 1e-5
#' @param maxiter the maximum number of iterations 
#'@return Returns a matrix of coefficients 
#'@author Brad Price <brad.price@@mail.wvu.edu>
#'@useDynLib mcen BinUp
bin_horse<-function(Y,X,delta,gamma_y,y_clusters,set_length,eps,maxiter){
  Beta01=matrix(0,nrow=dim(X)[2],ncol=dim(Y)[2])
  BetaI2=Beta01+1
  iter=0
  while(sum((Beta01-BetaI2)^2)/length(Beta01)>eps && iter<maxiter){
    BetaI2=Beta01
    Mine<-CalcHorseBin(Y,X,BetaI2)
    Out<-.C("BinUp",TmV=as.double((as.vector(Mine$TmV))),Tm=as.double(as.vector(Mine$Tm)),V=as.double(as.vector(Mine$V)),C=as.double(as.vector(Mine$C)),R=as.double(as.vector(Mine$R)),
            Beta=as.double(as.vector(BetaI2)),r=as.integer(dim(Y)[2]),p=as.integer(dim(X)[2]),n=as.integer(dim(X)[1]),
            y_clusters=as.double(as.vector(y_clusters)), sl=as.double(as.vector(set_length)), gamma_y=as.double(gamma_y), delta=as.double(delta),
            iter=as.double(maxiter),eps=as.double(eps),third=as.double(as.vector(Beta01)),fourth=as.double(as.vector(Beta01)),
            sixth=as.double(as.vector(Beta01)),first=as.double(as.vector(Beta01)),second=as.double(as.vector(Beta01)),fifth=as.double(as.vector(Beta01)))
    Beta01=matrix(Out$Beta,byrow = FALSE,ncol=dim(Y)[2],nrow=dim(X)[2])
    ##print(sum((Beta01-BetaI2)^2)/length(Beta01))
    iter=iter+1
  }
  ##dyn.unload("BinUpdate.so")
  return(list(Beta=Beta01))
}


#' Calculates cluster assignment and coefficient estimates for a binomial mcen.
#'@param beta Initial estimate of coefficients.
#'@param delta Tuning parameter for L1 penalty.
#'@param y Matrix of responses.
#'@param x Matrix of predictors. 
#'@param family type of likelihood used two options "mgaussian" or "mbinomial"
#'@param ky Number of clusters used for grouping response variables. 
#'@param gamma_y Tuning parameter for the penalty between fitted values for responses in the same group.
#'@param eps Convergence criteria
#'@param clusterMethod Which clustering method was used, currently support kmeans or kmeanspp
#'@param clusterIterations Number of iterations for cluster convergence
#'@param clusterStartNum Number of random starting points for clustering
#'@param cluster_y An a priori definition of clusters. If clusters are provided they will remain fixed and are not estimated. Objective function is then convex. 
#'@param max_iter The maximum number of iterations for estimating the coefficients
#'@author Brad Price <brad.price@@mail.wvu.edu> 
mcen_bin_workhorse<-function(beta,delta=NULL,y,x,family="mbinomial",ky=NULL,gamma_y=1,eps=.00001,
                         clusterMethod="kmeans",clusterIterations=100,clusterStartNum=30,cluster_y=NULL,max_iter=10){
  beta_dim <- dim(beta)
  p <- beta_dim[1]
  r <- beta_dim[2]
  if(sum(beta==0)==p*r){
    #if initial estimate is totally sparse then mcen update is not needed. Philosophically, I'm ok with this because mcen idea is 
    # to group coefficient estimates, if all coefficients are zero there is nothing to group and nothing for mcen to do. 
    return_val <- list(beta = as(beta,"sparseMatrix"), clusters=rep(0,r))
    return( return_val)
  } else{
    fitted_vals <- x%*%beta
    fitted_vals<-as.matrix(fitted_vals)
    iter <- 0
    fitting_done <- FALSE
    while(!fitting_done){
      #print(paste("fitting done is ", fitting_done, " and iter is ", iter))
      num_potential_centers <- ncol(unique(fitted_vals,MARGIN = 2))
      if(num_potential_centers < ky){
        fitting_done <- TRUE
        warning(paste0("Algorithm stopped early for gamma_y ", gamma_y, ", delta ", delta, " and ky ", ky, " because number of unqiue fitted values is less than ", ky, " Common for large values of deltas that cause spares solution and little variability in fitted values.\n")) 
        clusters_used <- rep(0,r)
      } else if(is.null(cluster_y)==FALSE){
        clusters_used <- cluster_y
      } else{
        y_clusters <- cluster(fitted_vals,ky,clusterMethod,clusterIterations,clusterStartNum)
        clusters_used <- y_clusters$cluster
      }
      if(iter > 0){
        if(is.null(cluster_y) && fitting_done==FALSE){
          #fitting_done=SetEq(y_clusters2$centers,y_clusters$centers)
          fitting_done=SetEq(y_clusters2$cluster,y_clusters$cluster)
        } else{
          fitting_done <- TRUE
        }
      }
      if(fitting_done == FALSE){
      sl=NULL
      for(m in 1:length(clusters_used)){
        sl[m]=length(which(clusters_used==clusters_used[m]))
      }
      beta<-bin_horse(Y=y,X=x,delta = delta,gamma_y = gamma_y,y_clusters =clusters_used ,set_length = sl,eps=eps,maxiter = max_iter)$Beta
      } 
      if(iter == max_iter){
        fitting_done <- TRUE
      }
      iter <- iter + 1
      if(is.null(cluster_y)==TRUE && fitting_done==FALSE){
        y_clusters2 <- y_clusters
      }
      fitted_vals=x%*%beta
    }
    if(is.null(cluster_y)==FALSE){
      fitting_done <- TRUE #only one iteration of fitting done if clusters are prespecified
    }
  }
  return_val <- list(beta = as(beta,"sparseMatrix"), clusters=clusters_used)
  return( return_val)
}							


