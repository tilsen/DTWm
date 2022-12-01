
get_pred <- function(model,times,conds,subjs) {
  #generate predictions for export
  
  nconds = length(conds)
  
  newd <- expand.grid(0,conds,times)
  newd <- data.frame(subj=newd[[1]],cond=newd[[2]],time=newd[[3]])
  
  #note: returns 1.96*s.e. 
  p <- get_predictions(model,
                       cond=list(subj=subjs[1],cond=conds,time=times),
                       se=TRUE,
                       rm.ranef = TRUE,print.summary = FALSE)
  newd$mu <- p$fit
  newd$ci <- p$CI
  
  for (i in (1:length(subjs))){
    newd$subj <- i

    #note: returns s.e.
    p <- predict.gam(model,newd,
                     se.fit=TRUE,
                     print.summary = FALSE)
    varname_mu <- paste0('mu_s',i)
    varname_ci <- paste0('ci_s',i)
    newd[[varname_mu]] <- p$fit
    newd[[varname_ci]] <- p$se.fit  
  }
  
  newdd <- data.frame(time=times)
  
  #difference predictions across subj
  for (i in 1:(nconds-1)) {
    for (j in (i+1):nconds){

      #returns 1.96*s.e.
      dd <- get_difference(model,
                           comp=list(cond=c(i,j)),
                           cond=list(time=times),
                           rm.ranef=TRUE,
                           se=TRUE,print.summary = FALSE)
      varname_diff <- paste0("d_",i,"_",j)
      varname_ci <- paste0("ci_",i,"_",j)
      newdd[[varname_diff]] <- dd$difference
      newdd[[varname_ci]] <- dd$CI
    }
  }
  
  #difference predictions by subj
  for (s in (1:length(subjs))){
    for (i in 1:(nconds-1)) {
      for (j in (i+1):nconds){

        #returns 1.96*s.e.
        dd <- get_difference(model,
                             comp=list(cond=c(i,j)),
                             cond=list(subj=s,time=times),
                             rm.ranef=FALSE,
                             se=TRUE,print.summary = FALSE)
        varname_diff <- paste0("d_",i,"_",j,"_",s)
        varname_ci <- paste0("ci_",i,"_",j,"_",s)
        newdd[[varname_diff]] <- dd$difference
        newdd[[varname_ci]] <- dd$CI
      }
    }     
  }
  
  pred <- rbind.fill(newd,newdd)
  
  return(pred)
}