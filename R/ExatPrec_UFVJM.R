#'@export
#'@description
#'#Algortimo para avaliar a exatidao e precisao dos dados
#' @param x Concentracao a ser avaliada
#' @param y respostas cromatográficas de x
#' @param z é a resposta cromatográfica ...
#'@example
#'x<- 0.02
#'y<- c(233153.0,257269.0,307816.0,256882.0,291418.0,264616.0,262520.0)
#'z<- 267896
#'EP(x,y,z)


# Algoritmo para verificar a exatidao e precisao
EP<- function(x,y,z){
  dpr<- (sd(y)/median(y)*100)
  rec<- (median(y)/z) *100
  cat("Resultado:",rec,"+-",dpr)

  cat("\nSe seu resultado estiver entre 70 e 120 & seu desvio for menor que 20.
Seu resultado é preciso e exato")
}
