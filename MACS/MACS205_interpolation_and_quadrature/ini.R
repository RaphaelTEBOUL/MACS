##Création des histogramme basés sur la loi beta et la loi de X
rm(list=ls())
graphics.off()

##Nombre d'échantillons

n_sample <- 10^(7)

generation_X = function(n_sample)
{
  X = rbeta(n = n_sample, apar, bpar)
  return(X)
}


hist_value = function(histog, x)
{
  if(any(x<0)|| any( x>=1))
  {
    stop('x n est pas dans [0,1]')
  }
  else{
    index_x <- floor(x/space_int);
    y <- histog[index_x+1];
    return(y)
  }
}



## parametres globaux

space_int <- 2^(-8)
## largeur des barres de l'histogramme 

bin_vec <- seq(space_int/2, (1 - space_int/2), by=space_int)
## vecteur des centres des barres de l'histogrammes

apar <- 1.7 ## parametre alpha de la loi beta
bpar <- 5.1 ## parametre alpha de la loi beta



##Génération de  l'histogramme
X_sample = generation_X(n_sample);
hist1 <- hist(X_sample, breaks = seq(0,1,by=space_int), probability= FALSE)$counts/ (space_int*n_sample);
## hist1 est un vecteur de même taille que bin_vec . 

rm(X_sample)

##Tracé de l'histogramme


plot(bin_vec,hist1, type='h') 
save(hist1, file="hist1.rda")
##source("functions_ETC2018.R")



