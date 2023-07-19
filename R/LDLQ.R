#' @title LDLQ
#' @name Limite de detecção e quantificação
#'
#' @description A função informa o Limite de Detecção do analito e seu Limite de
#' Quantificação
#'
#' @param x Valor da área fortificada
#' @param y valor da área do branco
#'
#' @return Valores de LD e LQ
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#'x = 49.480
#'y = 17.054
#'LDLQ(x,y)
#'
#'@export
#Função LDLQ simplificada
LDLQ<- function(x,y){
  LD = (x/y) * 3 #x é a área fortificada e y é a área do branco
  LQ = (x/y) * 10
  cat("Limite de detecção:",LD)
  cat("\nLimite de quantificação:",LQ)
}
