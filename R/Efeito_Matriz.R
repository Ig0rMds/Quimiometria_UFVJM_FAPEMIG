#' @title EM
#' @name Efeito de Matriz
#'
#' @description A função analisa o efeito de matriz presente da comparação da amostra
#' obtida com a análise feita em branco
#'
#' @param x Um vetor da variável independente X
#' @param y Um vetor da variável dependente Y, sendo a resposta da amostra fortificada
#' @param z Um vetor da variável dependente Z, sendo a resposta do branco
#'
#' @return Avaliação a respeito do efeito de matriz em porcentagem
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#' #teste com TCDF pelo MELDD (Dados obtidos de SICUPIRA,2019)
#'x<- c(5.3,5.3,5.3,13.3,13.3,13.3,21.3,21.3,21.3,29.3,29.3,29.3,37.3,37.3,37.3,45.3,45.3,45.3)
#'y<- c(423832,267896,370213,1545854,1093452,780549,2378730,2204902,2709421,3492212,3591964,3523445,4762476,4830668,4853425,5334327,5753813,5648307)
#'z<- c(733485,755427,760125,1827101,1860320,1853241,2904764,2913248,3025938,4418976,4492318,4516138,5722249,5661746,5535078,6741298,7079091,7005771)
#'EM(x,y,z)
#'
#'
#'@export
#Pacote Quimiometria_UFVJM_FAPEMIG
#Algortimo para Efeito de Matriz
###################################################################
EM<- function (x,y,z){
  library(car)
  #fazendo a regressão
  regmat = lm(y~x) #regressao da matriz
  regsol = lm(z~x) #regressao do solvente
  #coeficientes lineares
  b0mat<-coefficients(regmat)[2]
  b0sol<-coefficients(regsol)[2]
  #Calculando o efeito de matriz
  Res <- (b0mat/b0sol)*100
  #condicoes
  if (Res < 100) {
    cat("Resultado:",Res,"% Logo, há uma reducao na resposta cromatográfica\n------------------------------------------------------------------------\n")
  }
  if (Res > 100) {
    cat("Resultado:",Res,"% Logo, há um aumento na resposta cromatográfica\n------------------------------------------------------------------------\n")
  }
  if (Res == 100) {
    cat("Resultado:",Res,"% Logo, Não há efeito de matriz\n------------------------------------------------------------------------\n")
  }
}
