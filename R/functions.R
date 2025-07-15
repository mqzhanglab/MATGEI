#' run multinomial regression and calc the pvalues of invariant
#'
#' this function  calculate the core pvalues of invariant. This function is wraped by run_invariant_pval
#' This function construct the input for function run_invariant_pval
#' @param df_pheno the matrix that each row is a sample and each column is a feature.
#' @param df_geno the matrix of geno information, the additive model with each column is a sample, each row is a SNP.
#' @param df_pc10 the matrix of first 10 pcs from geno, in the 10 columns.
#' @param X_tem column names in df_pheno that used for all covariates(non-interest items)
#' @param G_term column names in df_pheno that used for all Geno parts GxE parts(interest items to be test)
#' @param Y_term column names in df_pheno that used for interested response.
#'
#' @return
#' @export
#'
#' @examples
#' X_item=c("sex","age",
#' "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
#' E_item=c("sun_hour_summer")
#' G_item=c("geno")
#' GE_item=c("sun_hour_summer:geno")
#' Y_item=c("response")
#' run_invariant_pval(df_pheno=cur_pheno,df_geno=geno,df_pc10=pc10_geno,
#'                X_term=c(X_item,E_item),G_term=c(G_item,GE_item),Y_term=Y_item)
#'
run_invariant_pval=function(df_pheno,df_geno,df_pc10=NA,X_term,G_term,Y_term){
  all_dat=cbind(df_pheno,df_pc10)

  #the computation time are similar if we compare the following two options and only have one core.
  #cal option 1
  reslist=list()
  for(i in 1:nrow(df_geno)){
    gc()
    all_dat$geno=df_geno[i,]
    #exclude SNPs that >10% are missing or non-missing record <10
    if(sum(!is.na(all_dat$geno))<10 | sum(is.na(all_dat$geno))/length(all_dat$geno)>0.1 ){
      reslist[[i]]=NA
    }else{
      cur_dat=all_dat[!is.na(all_dat$geno),]
      X=model.matrix(as.formula(paste("~", paste(X_term,collapse= "+"))),data=cur_dat)
      Y=model.matrix(as.formula(paste("~0+", paste(Y_term,collapse= "+"))),data=cur_dat)
      G=model.matrix(as.formula(paste("~0+", paste(G_term,collapse= "+"))), data=cur_dat)

      cur_res=tryCatch({nnet::multinom(as.formula(paste(paste(Y_term,collapse= "+"),"~", paste(X_term,collapse= "+"))),data=cur_dat)}, error = function(e) {NA})
      mu=tryCatch({fitted(cur_res)}, error = function(e) {NA})
      reslist[[i]]= tryCatch({invariant_pval_cal2(X,G,Y,mu)}, error = function(e) {rep(NA,8)}) #!! no tryCatch here.
    }
    print(paste("current SNP: ",i,Sys.time()))

    #print(reslist[[i]])
  }

  # #cal option 2
  # `%dorng%` <- doRNG::`%dorng%`
  # `%dopar%` <- foreach::`%dopar%`
  # reslist=foreach::foreach(i=1:nrow(df_geno)) %dorng% {
  #   #for(ip in 0:10) {
  #   #print(ip)
  #   all_dat$geno=df_geno[i,]
  #   cur_dat=all_dat[!is.na(all_dat$geno),]
  #   X=model.matrix(as.formula(paste("~", paste(X_term,collapse= "+"))),data=cur_dat)
  #   G=model.matrix(as.formula(paste("~0+", paste(G_term,collapse= "+"))), data=cur_dat)
  #   Y=model.matrix(as.formula(paste("~0+", paste(Y_term,collapse= "+"))),data=cur_dat)
  #   cur_res=nnet::multinom(as.formula(paste(paste(Y_term,collapse= "+"),"~", paste(X_term,collapse= "+"))),data=cur_dat)
  #   mu=fitted(cur_res)
  #   invariant_pval_cal(X,G,Y,mu)

  names(reslist)=rownames(df_geno)
  return(reslist)
}

#' run multinomial regression and calc the pvalues of invariant for G+GT
#'
#' The difference between this function with run_invariant_pval is this function only calc F one time.
#' this function  calculate the core pvalues of invariant. This function is wraped by run_invariant_pval
#' This function construct the input for function run_invariant_pval
#' @param df_pheno the matrix that each row is a sample and each column is a feature.
#' @param df_geno the matrix of geno information, the additive model with each column is a sample, each row is a SNP.
#' @param df_pc10 the matrix of first 10 pcs from geno, in the 10 columns.
#' @param X_tem column names in df_pheno that used for all covariates(non-interest items)
#' @param G_term column names in df_pheno that used for all Geno parts GxE parts(interest items to be test)
#' @param Y_term column names in df_pheno that used for interested response.
#'
#' @return
#' @export
#'
#' @examples
#' X_item=c("sex","age",
#' "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
#' E_item=c("sun_hour_summer")
#' G_item=c("geno")
#' GE_item=c("sun_hour_summer:geno")
#' Y_item=c("response")
#' run_invariant_pval(df_pheno=cur_pheno,df_geno=geno,df_pc10=pc10_geno,
#'                X_term=c(X_item,E_item),G_term=c(G_item,GE_item),Y_term=Y_item)
#'
run_invariant_pval_GGT=function(df_pheno,df_geno,df_pc10=NA,X_term,G_term,Y_term){
  all_dat=cbind(df_pheno,df_pc10)

  #the computation time are similar if we compare the following two options and only have one core.
  #cal option 1
  reslist=list()

  #calculate the background info and F
  X=model.matrix(as.formula(paste("~", paste(X_term,collapse= "+"))),data=all_dat)
  Y=model.matrix(as.formula(paste("~0+", paste(Y_term,collapse= "+"))),data=all_dat)
  cur_res=tryCatch({nnet::multinom(as.formula(paste(paste(Y_term,collapse= "+"),"~", paste(X_term,collapse= "+"))),data=all_dat)}, error = function(e) {NA})
  mu_base=tryCatch({fitted(cur_res)}, error = function(e) {NA})

  mu_base=t(mu_base)
  J=nrow(mu_base)
  n=ncol(mu_base)

  #step 1: get V
  #print("cal Fu_list")
  Fu_list=array(dim=c(n*J,n*J))

  for(l in 1:J){
    for(t in 1:l){
      if(l==t){
        cur_Fu=diag(n)*(mu_base[l,]*(1-mu_base[l,]))
      }
      if(l!=t){
        cur_Fu=diag(n)*(-mu_base[l,]*mu_base[t,])
      }
      Fu_list[seq((l-1)*n+1,(l-1)*n+n),seq((t-1)*n+1,(t-1)*n+n)]=cur_Fu
      Fu_list[seq((t-1)*n+1,(t-1)*n+n),seq((l-1)*n+1,(l-1)*n+n)]=cur_Fu
    }
  }
  pval_sep=list()
  pval_sep_liu=list()
  pval_cauchy=list()
  pval_cauchy_liu=list()
  pval_integ=list()
  pval_integ_liu=list()
  V_sep=list()
  Q_sep=list()
  for(i_j in 1:J){
    gc()

    X_bar= kronecker(diag(J-1), X) #kronecker product
    FuJ=Fu_list[-seq((i_j-1)*n+1,(i_j-1)*n+n),-seq((i_j-1)*n+1,(i_j-1)*n+n)]
    #FuJ_solve1=tryCatch({solve(t(X_bar) %*% FuJ %*% X_bar )}, error = function(e) {NA})
    ## V_adjust=FALSE
    #if(is.na(Fu_solve1)){
      FuJ_solve1=solve(t(X_bar) %*% FuJ %*% X_bar + (1e-16)*diag(rep(1, ncol(X_bar))))
    #  # V_adjust=TRUE
    #}
    FuJ_solve2=FuJ %*% X_bar %*% FuJ_solve1 %*% t(X_bar)%*% FuJ
    FuJ_ready=(FuJ-FuJ_solve2)

    #add on each gene's situations.
    for(i in 1:nrow(df_geno)){
      gc()
      if(i_j==1){
        pval_sep[[i]]=list()
        pval_cauchy[[i]]=matrix(ncol=1,nrow=J)
        pval_cauchy_liu[[i]]=matrix(ncol=1,nrow=J)
        V_sep[[i]]=list()
        Q_sep[[i]]=list()
      }
      all_dat$geno=df_geno[i,]

      #exclude SNPs that >10% are missing or non-missing record <10
      if(sum(!is.na(all_dat$geno))<10 | sum(is.na(all_dat$geno))/length(all_dat$geno)>0.1 ){
        reslist[[i]]=NA
      }else{
        cur_index=which(!is.na(all_dat$geno))
        cur_dat=all_dat[cur_index,]
        mu=mu_base[,cur_index] # we will calculate mu based on a larger matrix
        X=model.matrix(as.formula(paste("~", paste(X_term,collapse= "+"))),data=cur_dat)
        Y=model.matrix(as.formula(paste("~0+", paste(Y_term,collapse= "+"))),data=cur_dat)
        G=model.matrix(as.formula(paste("~0+", paste(G_term,collapse= "+"))), data=cur_dat)
        X_bar= kronecker(diag(J-1), X) #kronecker product
        G_bar= kronecker(diag(J-1), G)

        p=ncol(G)

        keep_index=rep(cur_index,times=(J-1))+n*rep(0:(J-2),each=length(cur_index))
        if(n==length(cur_index)){
          Fu_ready=FuJ_ready
        }else{
          Fu=FuJ[keep_index,keep_index]
          # #option2: get V from designed method.
          #Fu_solve1=tryCatch({solve(t(X_bar) %*% Fu %*% X_bar )}, error = function(e) {NA})
          #V_adjust=FALSE
          #if(is.na(Fu_solve1)){
            Fu_solve1=solve(t(X_bar) %*% Fu %*% X_bar + (1e-16)*diag(rep(1, ncol(X_bar))))
          #  V_adjust=TRUE
          #}
          Fu_solve2=Fu %*% X_bar %*% Fu_solve1 %*% t(X_bar)%*% Fu
          Fu_ready=(Fu-Fu_solve2)
        }
        #Method1_cauchyP

        cur_V=t(G_bar) %*% Fu_ready %*% G_bar
        # ## ###
        #Step2: Get Q
        S=t(G) %*% (Y-t(mu))

        cur_Q=sum(S[,-i_j]^2)
        V_sep[[i]][[i_j]]=cur_V
        Q_sep[[i]][[i_j]]=cur_Q
        #run davies method
        pval_sep[[i]][[i_j]]=tryCatch({SKAT:::Get_PValue(K=cur_V,Q=cur_Q)}, error = function(e) {NA})
        pval_cauchy[[i]][i_j]=tryCatch({pval_sep[[i]][[i_j]]$p.value}, error = function(e) {NA})
        pval_cauchy_liu[[i]][i_j]=tryCatch({pval_sep[[i]][[i_j]]$p.val.liu}, error = function(e) {NA})

        #Method2_integrative
        # just take the last V
        if(i_j==J){
          D=kronecker(rbind(diag(J-1),-1),diag(p))
          cur_L=D %*% cur_V %*% t(D)
          pval_combine=SKAT:::Get_PValue(K=cur_L,Q=sum(S^2))

          pval_integ[i]=tryCatch({pval_combine["p.value"]}, error = function(e) {NA})
          pval_integ_liu[i]=tryCatch({pval_combine["p.val.liu"]}, error = function(e) {NA})
        }
      }
      print(paste("current category: ",i_j,", current SNP: ",i,Sys.time()))
    }
  }
  for(i in 1:nrow(df_geno)){
    pval_cauchy[[i]]=tryCatch({ACAT::ACAT(pval_cauchy[[i]])}, error = function(e) {NA})
    pval_cauchy_liu[[i]]=tryCatch({ACAT::ACAT(pval_cauchy_liu[[i]])}, error = function(e) {NA})

    #print(c("pval_sep 1: ",pval_sep[[1]]))
    reslist[[i]]=list(pval_cauchy=pval_cauchy[[i]],
                      pval_cauchy_liu=pval_cauchy_liu[[i]],
                      pval_integ=pval_integ[i],
                      pval_integ_liu=pval_integ_liu[i],
                      pval_sep=sapply(pval_sep[[i]],function(x)x["p.value"]),
                      pval_sep_liu=sapply(pval_sep[[i]],function(x)x["p.val.liu"]),
                      V_sep=V_sep[[i]],
                      Q_sep=Q_sep[[i]])
  }
  names(reslist)=rownames(df_geno)
  return(reslist)
}

#' this function calculate the core pvalues of invariant. This function is wraped by run_invariant_pval
#'
#' @param X matrix contains all covariates(non-interest items)
#' @param G matrix contains all Geno parts GxE parts(interest items to be test)
#' @param Y matrix contains interested response.
#' @param mu estimated y from  multinomial regression. fitted values of output of nnet::multinom
#' @param VQstat bool, if TRUE, returns the V & Q statistics as well.
#'
#' @return
#' @export
#'
#' @examples
invariant_pval_cal=function(X,G,Y,mu,VQstat=FALSE){
  J=ncol(mu)
  n=nrow(mu)
  p=ncol(G)
  #step 1: get V
  #print("cal Fu_list")
  Fu_list=array(dim=c(n*J,n*J))

  for(l in 1:J){
    for(t in 1:l){
      if(l==t){
        cur_Fu=diag(mu[,l]*(1-mu[,l]))
      }
      if(l!=t){
        cur_Fu=diag(-mu[,l]*mu[,t])
      }
      Fu_list[seq((l-1)*n+1,(l-1)*n+n),seq((t-1)*n+1,(t-1)*n+n)]=cur_Fu
      Fu_list[seq((t-1)*n+1,(t-1)*n+n),seq((l-1)*n+1,(l-1)*n+n)]=cur_Fu
    }
  }

  #Step2: Get Q
  S=t(G) %*% (Y-mu)

  pval_sep=list()
  pval_cauchy=matrix(ncol=1,nrow=J)
  pval_cauchy_liu=matrix(ncol=1,nrow=J)
  V_sep=list()
  Q_sep=list()
  #Method1_cauchyP
  for(i_j in 1:J){
    gc()
    #cur_V=get_V(X,G,Fu_list,i_rm=i_j,all_cat=FALSE)
    ## ####
    X_bar= kronecker(diag(J-1), X) #kronecker product
    G_bar= kronecker(diag(J-1), G)
    Fu=Fu_list[-seq((i_j-1)*n+1,(i_j-1)*n+n),-seq((i_j-1)*n+1,(i_j-1)*n+n)]


    # #print(c("cal cur_V "))
    # #option 1: get V from eigen for computational benefits.
    # eigenFu=eigen(-Fu)
    # min(eigenFu$values)
    # eigenFu$values=pmax(0,eigenFu$values)
    # L=sqrt(eigenFu$values) * t(eigenFu$vector)
    # M=L%*%X_bar
    # N=L%*%G_bar
    #
    # ## pivoted Chol for `M'M`; we want lower triangular factor `L = R'`:
    # ## we also suppress possible rank-deficient warnings (no harm at all!)
    # LL=t(suppressWarnings(chol(crossprod(M), pivot = TRUE)))
    # ## compute Qt
    # r = attr(LL, "rank")
    # piv = attr(LL, "pivot")
    # Qt = forwardsolve(LL, t(M[, piv]), r)
    # Hat = crossprod(Qt)
    # diag(Hat)=diag(Hat)-1
    # cur_V=t(N)%*% Hat %*% N

    # #option2: get V from designed method.
    #Fu_solve1=tryCatch({solve(t(X_bar) %*% Fu %*% X_bar )}, error = function(e) {NA})
    #V_adjust=FALSE
    #if(is.na(Fu_solve1)){
      Fu_solve1=solve(t(X_bar) %*% Fu %*% X_bar + (1e-16)*diag(rep(1, ncol(X_bar))))
    #  V_adjust=TRUE
    #}
    Fu_solve2=Fu %*% X_bar %*% Fu_solve1 %*% t(X_bar)%*% Fu
    cur_V=t(G_bar) %*% (Fu-Fu_solve2) %*% G_bar
    # ## ###

    cur_Q=sum(S[,-i_j]^2)
    V_sep[[i_j]]=cur_V
    Q_sep[[i_j]]=cur_Q
    #run davies method
    pval_sep[[i_j]]=tryCatch({SKAT:::Get_PValue(K=cur_V,Q=cur_Q)}, error = function(e) {NA})
    pval_cauchy[i_j]=tryCatch({pval_sep[[i_j]]$p.value}, error = function(e) {NA})
    pval_cauchy_liu[i_j]=tryCatch({pval_sep[[i_j]]$p.val.liu}, error = function(e) {NA})
  }

  pval_cauchy=tryCatch({ACAT::ACAT(pval_cauchy)}, error = function(e) {NA})
  pval_cauchy_liu=tryCatch({ACAT::ACAT(pval_cauchy_liu)}, error = function(e) {NA})

  #Method2_integrative
  # just take the last V
  D=kronecker(rbind(diag(J-1),-1),diag(p))
  cur_L=D %*% cur_V %*% t(D)
  pval_combine=SKAT:::Get_PValue(K=cur_L,Q=sum(S^2))

  pval_integ=tryCatch({pval_combine["p.value"]}, error = function(e) {NA})
  pval_integ_liu=tryCatch({pval_combine["p.val.liu"]}, error = function(e) {NA})

  #print(c("pval_sep 1: ",pval_sep[[1]]))
  if(VQstat==FALSE){
    res=list(pval_cauchy=pval_cauchy,
             pval_cauchy_liu=pval_cauchy_liu,
             pval_integ=pval_integ,
             pval_integ_liu=pval_integ_liu,
             pval_sep=sapply(pval_sep,function(x)x["p.value"]),
             pval_sep_liu=sapply(pval_sep,function(x)x["p.val.liu"]))
  }
  if(VQstat==TRUE){
    res=list(pval_cauchy=pval_cauchy,
             pval_cauchy_liu=pval_cauchy_liu,
             pval_integ=pval_integ,
             pval_integ_liu=pval_integ_liu,
             pval_sep=sapply(pval_sep,function(x)x["p.value"]),
             pval_sep_liu=sapply(pval_sep,function(x)x["p.val.liu"]),
             V_sep=V_sep,
             Q_sep=Q_sep)
  }
  #print(res)
  return(res)
}


#' this function calculate the core pvalues of invariant. This function is wraped by run_invariant_pval
#'
#' @param X matrix contains all covariates(non-interest items)
#' @param G matrix contains all Geno parts GxE parts(interest items to be test)
#' @param Y matrix contains interested response.
#' @param mu estimated y from  multinomial regression. fitted values of output of nnet::multinom
#' @param VQstat bool, if TRUE, returns the V & Q statistics as well.
#'
#' @return
#' @export
#'
#' @examples
invariant_pval_cal2=function(X,G,Y,mu,VQstat=FALSE){
  J=ncol(mu)
  n=nrow(mu)
  p=ncol(G)
  #step 1: get V
  # V_adjust=FALSE
  #V_sep=tryCatch({get_V3(X,G,Y,mu)}, error = function(e) {NA}) #RCPP version
  #if(is.na(V_sep)){
    V_sep=tryCatch({get_V3_1(X,G,Y,mu)}, error = function(e) {NA}) #RCPP version
    # V_adjust=TRUE
  #}

  #V_sep=get_V2(X,G,Y,mu) #R version
  #Step2: Get Q
  S=t(G) %*% (Y-mu)
  Q_sep=list()

  #step3:get pval
  pval_sep=list()
  pval_cauchy=matrix(ncol=1,nrow=J)
  pval_cauchy_liu=matrix(ncol=1,nrow=J)

  #Method1_cauchyP
  for(i_j in 1:J){
    # ## ###
    cur_V=V_sep[[i_j]]
    cur_Q=sum(S[,-i_j]^2)
    Q_sep[[i_j]]=cur_Q
    #run davies method
    pval_sep[[i_j]]=tryCatch({SKAT:::Get_PValue(K=cur_V,Q=cur_Q)}, error = function(e) {NA})
    pval_cauchy[i_j]=tryCatch({pval_sep[[i_j]]$p.value}, error = function(e) {NA})
    pval_cauchy_liu[i_j]=tryCatch({pval_sep[[i_j]]$p.val.liu}, error = function(e) {NA})
  }

  pval_cauchy=tryCatch({ACAT::ACAT(pval_cauchy)}, error = function(e) {NA})
  pval_cauchy_liu=tryCatch({ACAT::ACAT(pval_cauchy_liu)}, error = function(e) {NA})

  #Method2_integrative
  # just take the last V
  D=kronecker(rbind(diag(J-1),-1),diag(p))
  cur_L=D %*% cur_V %*% t(D)
  pval_combine=SKAT:::Get_PValue(K=cur_L,Q=sum(S^2))

  pval_integ=tryCatch({pval_combine["p.value"]}, error = function(e) {NA})
  pval_integ_liu=tryCatch({pval_combine["p.val.liu"]}, error = function(e) {NA})
  #print(c("pval_sep 1: ",pval_sep[[1]]))
  if(VQstat==FALSE){
    res=list(pval_cauchy=pval_cauchy,
             pval_cauchy_liu=pval_cauchy_liu,
             pval_integ=pval_integ,
             pval_integ_liu=pval_integ_liu,
             pval_sep=sapply(pval_sep,function(x)x["p.value"]),
             pval_sep_liu=sapply(pval_sep,function(x)x["p.val.liu"]#,
             #V_adjust=V_adjust
             ))
  }
  if(VQstat==TRUE){
    res=list(pval_cauchy=pval_cauchy,
             pval_cauchy_liu=pval_cauchy_liu,
             pval_integ=pval_integ,
             pval_integ_liu=pval_integ_liu,
             pval_sep=sapply(pval_sep,function(x)x["p.value"]),
             pval_sep_liu=sapply(pval_sep,function(x)x["p.val.liu"]),
             #V_adjust=V_adjust,
             V_sep=V_sep,
             Q_sep=Q_sep)
  }



  #print(res)
  return(res)
}

#' get statistic V
#'
#' this function calculate the statistics V for each categories. This function is wraped by invariant_pval_cal
#' In practice, we usually use the Rcpp function get_V for computation speed benefit.
#'
#' @param X matrix contains all covariates(non-interest items)
#' @param G matrix contains all Geno parts GxE parts(interest items to be test)
#' @param Y matrix contains interested response.
#' @param mu estimated y from  multinomial regression. fitted values of output of nnet::multinom
#'
#' @return V_list a list of V's calculated.
#' @export
#'
#' @examples
get_V2=function(X,G,Y,mu){

  J=ncol(mu)
  n=nrow(mu)
  #step 1: get V
  #print("cal Fu_matrix")
  Fu_matrix=matrix(nrow=n*J,ncol=n*J)

  for(l in 1:J){
    for(t in 1:l){
      if(l==t){
        Fu_matrix[seq((l-1)*n+1,(l-1)*n+n),seq((t-1)*n+1,(t-1)*n+n)]=diag(mu[,l]*(1-mu[,l]))
        Fu_matrix[seq((t-1)*n+1,(t-1)*n+n),seq((l-1)*n+1,(l-1)*n+n)]=diag(mu[,l]*(1-mu[,l]))
      }
      if(l!=t){
        Fu_matrix[seq((l-1)*n+1,(l-1)*n+n),seq((t-1)*n+1,(t-1)*n+n)]=diag(n)*(-mu[,l]*mu[,t])
        Fu_matrix[seq((t-1)*n+1,(t-1)*n+n),seq((l-1)*n+1,(l-1)*n+n)]=diag(n)*(-mu[,l]*mu[,t])
      }
    }
  }

  V_list=list()
  #Method1_cauchyP
  for(i_j in 1:J){
    #cur_V=get_V(X,G,Fu_matrix,i_rm=i_j,all_cat=FALSE)
    ## ####
    X_bar= kronecker(diag(J-1), X) #kronecker product
    G_bar= kronecker(diag(J-1), G)
    if(i_j==1){
      Fu=Fu_matrix[seq((n+1),J*n),seq((n+1),J*n)]
    }else if(i_j==J){
      Fu=Fu_matrix[seq(1,(J-1)*n),seq(1,(J-1)*n)]
    }else if(i_j!=1 && i_j!=J ){
      Fu=Fu_matrix[c(seq(1,(i_j-1)*n),seq((i_j)*n+1,J*n)),c(seq(1,(i_j-1)*n),seq((i_j)*n+1,J*n))]
    }

    # #option1: get V from hat matrix.
    #L = tryCatch({chol(Fu)}, error = function(e) {NA})
    ## V_adjust=FALSE
    #if(is.na(L)){
      L = chol(Fu + (1e-16)*diag(rep(1, ncol(Fu))))
    #  #V_adjust=TRUE
    #}
    M=L%*%X_bar
    N=L%*%G_bar
    ## pivoted Chol for `M'M`; we want lower triangular factor `L = R'`:
    ## we also suppress possible rank-deficient warnings (no harm at all!)
    LL=t(suppressWarnings(chol(crossprod(M), pivot = TRUE)))
    ## compute Qt
    r = attr(LL, "rank")
    piv = attr(LL, "pivot")
    Qt = forwardsolve(LL, t(M[, piv]), r)
    Hat = crossprod(Qt)
    diag(Hat)=diag(Hat)-1
    cur_V=t(N)%*% (-Hat) %*% N

    # #option2: get V from designed method.
    # Fu_solve1=tryCatch({solve(t(X_bar) %*% Fu %*% X_bar )}, error = function(e) {NA})
    # V_adjust=FALSE
    # if(is.na(Fu_solve1)){
    #   Fu_solve1=solve(t(X_bar) %*% Fu %*% X_bar + (1e-16)*diag(rep(1, ncol(X_bar))))
    #   V_adjust=TRUE
    # }
    # Fu_solve2=Fu %*% X_bar %*% Fu_solve1 %*% t(X_bar)%*% Fu
    # cur_V=t(G_bar) %*% (Fu-Fu_solve2) %*% G_bar
    # ## ###

    V_list[[i_j]]=cur_V
  }
  return(V_list)
}

#' get the top PCs from PCA, Part I
#'
#' this function get the top PCs from PCA
#' @param count_m a matrix/data frame with each column as genes, each rows as subjects.
#' @param subj_id a vector with length as same as rows of count_m, it gives the id of subject of count_m
#' @param num_top_PCs the number of top PCs that would be picked.
#'
#' @return pca_res output list of function pca_analysis
#' @export
#'
#' @examples
#'
pca_analysis = function(count_m, subj_id, num_top_PCs=10){
  pca_model = prcomp(count_m)
  #print(summary(pca_model))
  #print(fviz_eig(pca_model))

  # variance explained by each principal component
  VE = pca_model$sdev^2

  # proportion of variance explained by each principal component
  PVE = VE / sum(VE)
  #num_top_PCs = which(cumsum(PVE) >= 0.9) %>% min
  #num_top_PCs = 10
  cat("num_top_PCs =",num_top_PCs,"\n")
  cum_var = cumsum(PVE)[num_top_PCs]
  cat("Variance explained by top PCs =", cum_var, "\n")
  # extract top PCs data
  topPCs = (pca_model$x) %>%
    as.data.frame %>%
    select(all_of(paste0("PC", 1:num_top_PCs))) %>%
    mutate(USUBJID = subj_id) %>%
    select(USUBJID, everything())

  return(list("topPCs"=topPCs, "pca_model"=pca_model, "num_top_PCs"=num_top_PCs, "cum_var"=cum_var))
}

#' get the top PCs from PCA ,Part II
#'
#' this function get the top PCs from PCA
#' @param pca_res output list of function pca_analysis
#' @param subj_id a vector with length as same as rows of count_m, it gives the id of subject of count_m
#' @param count_m  matrix/data frame with each column as genes, each rows as subjects.
#'
#' @return topPCs matrix with columns as expected columns of PCs.
#' @export
#'
#' @examples
predict_topPCs = function(pca_res, subj_id,count_m){
  pca_model = pca_res$pca_model
  num_top_PCs = pca_res$num_top_PCs
  pred_PCs = predict(pca_model, count_m)
  # extract top PCs data
  topPCs = (pred_PCs) %>%
    as.data.frame %>%
    select(all_of(paste0("PC", 1:num_top_PCs))) %>%
    mutate(USUBJID = subj_id) %>%
    select(USUBJID, everything())
  return(topPCs)
}

#' explore the data by NNet
#'
#' this function calculate the zscores and calculate the two-sided pvals from the results of nnet.
#' @param geno: the additive geno matrix with each column as a sample and each row as a variant.
#' @param pheno: the phenotype dataframe that each row is a sample.
#'        The column include a response column named as response. It is multinomial.While other columns are covariates or variables.
#' @param  nnet_obj: object of output of nnet::multinom
#' @return two-sided pvalues from zscores (coefficients/standard.error) nnet:: multinom regression.
#' @export
#'
#' @examples
cal_nnet_pval=function(nnet_obj){

  z=summary(nnet_obj)$coefficients/summary(nnet_obj)$standard.error
  pval=(0.5-abs(0.5-pnorm(z)))*2
  pval
}

#' pval_plot_hist
#' function pval_plot_hist
#' plot hist plot and QQ plot for output lists of total_pvals.
#' @param pval_list The list with each items are the matrix of pvalues from function cal_nnet_pval
#' @param flie_name The filename
#' @param main_name The mainname on the figures.
#' @return p save the plots with designed names.
#' @export
#'
#' @examples
#' pval_list=total_pval_uv
#' main_name="hist_pval_uv"
#' pval_plot(pval_list=total_pval_uv,
#' main_name="Distributions of pvalues with environment factor uv")
pval_plot_hist=function(pval_list, main_name){
  ng=length(pval_list)
  nr=dim(pval_list[[1]])[[1]]
  nb=dim(pval_list[[1]])[[2]]
  cur_df=data.frame(pval=unlist(pval_list),
                    response=rep(dimnames(pval_list[[1]])[[1]],times=ng*nb),
                    beta=rep(rep(dimnames(pval_list[[1]])[[2]],each=nr),times=ng))

  #remove NAs
  cur_df=cur_df[!is.na(cur_df$pval),]
  cur_df$pval[cur_df$pval==0]=min(cur_df$pval[cur_df$pval>0])

  p=ggplot(data = cur_df, aes(x = pval)) +
    geom_histogram(stat="bin",bins=50,position="identity") +
    ggtitle(main_name)+
    facet_grid(response~beta , scales = "free") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          #strip.background = element_blank(),
          plot.margin = unit(c(1,1,5,5), "mm"))
  p
  #ggexport(ggarrange(plotlist = p,ncol = 1,nrow =1),
  #         filename=paste0(file_name,"_",chr,".pdf"),width=5*nb,height=5*nr)
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

#' pval_plot_hist
#' function pval_plot_QQ
#' plot hist plot and QQ plot for output lists of total_pvals.
#' @param pval_list The list with each items are the matrix of pvalues from function cal_nnet_pval
#' @param flie_name The filename
#' @param main_name The mainname on the figures.
#' @return p save the plots with designed names.
#' @export
#'
#' @examples
#' pval_list=total_pval_uv
#' main_name="hist_pval_uv"
#' pval_plot_QQ(pval_list=total_pval_uv,
#' main_name="Distributions of pvalues with environment factor uv")
#'
pval_plot_QQ=function(pval_list, main_name){
  ng=length(pval_list)
  nr=dim(pval_list[[1]])[[1]]
  nb=dim(pval_list[[1]])[[2]]
  cur_df=data.frame(pval=unlist(pval_list),
                    response=rep(dimnames(pval_list[[1]])[[1]],times=ng*nb),
                    beta=rep(rep(dimnames(pval_list[[1]])[[2]],each=nr),times=ng))
  #remove NAs
  cur_df=cur_df[!is.na(cur_df$pval),]
  cur_df$pval[cur_df$pval==0]=min(cur_df$pval[cur_df$pval>0])
  #add qq calculation

  p=
    ggplot(data = cur_df, aes(sample = pval)) +
    #geom_histogram(stat="bin",bins=50,position="identity") +
    stat_qq(distribution = stats::qunif)+
    #stat_qq_line()+
    geom_abline()+
    scale_x_continuous(trans=reverselog_trans(base=10),
                       labels=trans_format("identity", function(x) -x)) +
    scale_y_continuous(trans=reverselog_trans(base=10),
                       labels=trans_format("identity", function(x) -x)) +
    ggtitle(main_name) +
    facet_grid(response~beta )
  p
  #ggexport(ggarrange(plotlist = p,ncol = 1,nrow =1),
  #         filename=paste0(file_name,"_",chr,".pdf"),width=5*nb,height=5*nr)
}

#' run multinomial regression and calc the pvalues of glm
#' this function  calculate the pvalues with glm.
#' @param df_pheno the matrix that each row is a sample and each column is a feature.
#' @param df_geno the matrix of geno information, the additive model with each column is a sample, each row is a SNP.
#' @param df_pc10 the matrix of first 10 pcs from geno, in the 10 columns.
#' @param X_tem column names in df_pheno that used for all covariates(non-interest items)
#' @param G_term column names in df_pheno that used for all Geno parts GxE parts(interest items to be test)
#' @param Y_term column names in df_pheno that used for interested response.
#' @param Y_combine_index vector of index that used for which categories response would be collapse in one category for glm.
#'
#' @return
#' @export
#'
#' @examples
#' X_item=c("sex","age",
#' "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
#' E_item=c("sun_hour_summer")
#' G_item=c("geno")
#' GE_item=c("sun_hour_summer:geno")
#' Y_item=c("response")
#' #' Y_item=NA
#' run_glm_pval(df_pheno=cur_pheno,df_geno=geno,df_pc10=pc10_geno,
#'                X_term=c(X_item,E_item),G_term=c(G_item,GE_item),Y_term=Y_item)
run_glm_pval=function(df_pheno,df_geno,df_pc10=NA,X_term,G_term,Y_term,Y_combine_index=NA){
  all_dat=cbind(df_pheno,df_pc10)

  #the computation time are similar if we compare the following two options and only have one core.
  #cal option 1
  reslist=matrix(ncol=1,nrow=nrow(df_geno))
  for(i in 1:nrow(df_geno)){
    gc()
    all_dat$geno=df_geno[i,]
    #exclude SNPs that >10% are missing or non-missing record <10
    if(sum(!is.na(all_dat$geno))<10 | sum(is.na(all_dat$geno))/length(all_dat$geno)>0.1 ){
      reslist[i]=NA
    }else{
      cur_dat=all_dat[!is.na(all_dat$geno),]
      X=model.matrix(as.formula(paste("~", paste(X_term,collapse= "+"))),data=cur_dat)
      Y=model.matrix(as.formula(paste("~0+", paste(Y_term,collapse= "+"))),data=cur_dat)
      G=model.matrix(as.formula(paste("~0+", paste(G_term,collapse= "+"))), data=cur_dat)

      if(ncol(Y)>2){
        if(is.na(Y_combine_index)){
          Y_combine_index=seq(1,ceiling((ncol(Y))/2))
        }
        Yg=as.factor(rowSums(Y[,Y_combine_index]))
      }else{
        Yg=as.factor(Y[,1])
      }

      fullmodel=tryCatch({glm(as.formula(paste("Yg~", paste(X_term,collapse= "+"),"+",paste(G_term,collapse= "+"))),
                              data=cur_dat,family="binomial")}, error = function(e) {NA})
      nullmodel=tryCatch({glm(as.formula(paste("Yg~", paste(X_term,collapse= "+"))),
                              data=cur_dat,family="binomial")}, error = function(e) {NA})
      reslist[i]=tryCatch({anova(fullmodel,nullmodel,test="LRT")$"Pr(>Chi)"[2]}, error = function(e) {NA})

    }
    if(i%%1000==0){print(paste("current SNP: ",i,Sys.time()))} # 100 per min for 4k samples.

    #print(reslist[[i]])
  }

  names(reslist)=rownames(df_geno)
  return(reslist)
}

