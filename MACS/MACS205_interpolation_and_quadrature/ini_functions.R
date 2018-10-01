generation_X = function(n_sample)
{
    X = rbeta(n = n_sample, apar, bpar)
    return(X)
}


hist_value = function(histog, x)
{
    if(any(x<0)|| any( x>1))
    {
      stop('x n est pas dans [0,1]')
    }
    else{
        index_x <- floor(x/space_int);
        y <- histog[index_x+1];
        return(y)
    }
}

