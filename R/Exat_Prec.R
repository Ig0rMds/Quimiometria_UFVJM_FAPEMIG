#' @title EP
#' @name Exatidão & Precisão
#'
#' @description A função analisa a exatidão e precisão dos dados fornecidos
#'
#' @param x O valor da concentração analisada
#' @param y Um vetor da variável dependente Y
#' @param z O valor do parametro para a variável dependente
#'
#' @return Avaliação a respeito da exatidão e precisão
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#'Teste com TCDF metodologia MELDD
#'x<- 0.02
#'y<- c(233153.0,257269.0,307816.0,256882.0,291418.0,264616.0,262520.0)
#'z<- 267896
#'EP(x,y,z)
#'@export
#'Pacote Quimiometria_UFVJM_FAPEMIG
#'Algoritmo para verificar a exatidao e precisao
EP<- function(x,y,z){
  dpr<- (sd(y)/median(y)*100) #Cálculo do desv padrao relativo
  rec<- (median(y)/z) *100
  cat("Resultado:",rec,"+-",dpr)
}
