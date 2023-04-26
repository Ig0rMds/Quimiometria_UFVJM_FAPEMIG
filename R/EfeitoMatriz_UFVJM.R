#'@export
#'@description
#'#Algortimo para Efeito de Matriz
#' @param x são as concentrações
#' @param y são as respostas cromatográficas do extrato em matriz
#' @param z são as respostas cromatográficas do extrato em solvente
#'@example
#'x<- c(5.3,5.3,5.3,13.3,13.3,13.3,21.3,21.3,21.3,29.3,29.3,29.3,37.3,37.3,37.3,45.3,45.3,45.3)
#'y<- c(423832,267896,370213,1545854,1093452,780549,2378730,2204902,2709421,3492212,3591964,3523445,4762476,4830668,4853425,5334327,5753813,5648307)
#'z<- c(733485,755427,760125,1827101,1860320,1853241,2904764,2913248,3025938,4418976,4492318,4516138,5722249,5661746,5535078,6741298,7079091,7005771)
#'EM<-(x,y,z)

#entrada de dados para extrato de matriz e solvente
EM<- function (x,y,z){
  library(car)
  regmat = lm(y~x)
  regsol = lm(z~x)
  b0mat<-coefficients(regmat)[2]
  b0sol<-coefficients(regsol)[2]
  Res <- (b0mat/b0sol)*100
  #condicoes
  if (Res < 100) {
    cat("Resultado:",Res,"% Logo, há uma redução na resposta cromatográfica\n------------------------------------------------------------------------\n")
  }
  if (Res > 100) {
    cat("Resultado:",Res,"% Logo, há um aumento na resposta cromatográfica\n------------------------------------------------------------------------\n")
  }
  if (Res == 100) {
    cat("Resultado:",Res,"% Logo, Não há efeito de matriz\n------------------------------------------------------------------------\n")
  }
}
