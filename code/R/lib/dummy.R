corr_two_nominal_vars = function(v, dat){
    v=as.character(v)
    x = as.factor(as.character(dat[,v[1]]))
    y = as.factor(as.character(dat[,v[2]]))
    dd = data.frame(x=x,y=y)
    dd = dd[rowSums(is.na(dd))==0,]
    sm = chisq.test(dd$y, dd$x) # get a p-value
    r= assocstats(xtabs(~x+y,data=dd))$cramer # get a R^2 value for categorical vars
    pval = sm$p.value
    list(r,pval)
}
