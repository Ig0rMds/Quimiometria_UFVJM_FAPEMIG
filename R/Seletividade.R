#' @title Sel
#' @name Seletividade
#'
#' @description A função fornece o gráfico da curva da amostra fortificada, bem como
#' a de sua amostra branca, para fins de comparação visual a respeito da seletividade
#' da análise
#'
#' @param t Total data points
#' @param xm X axys multiplier
#' @param ym Y axys multiplier
#' @param y  resposta do analito fortificado
#' @param z  resposta da amostra em branco
#'
#' @return Gráfico para a análise da seletividade
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#' É difícil explicar como utilzar por aqui, confira minha tese para entender a
#' implmentação
#'
#'@export
#Pacote Quimiometria_UFVJM_FAPEMIG
#Algoritmo para seletividade
Sel<- function(t,xm,ym,y,z){
  x <- c(1: t-1)
  xt<- xm*x
  yr<- y[1:t]
  zr<- z[1:t]
  yt<- ym*y
  zt<- ym*z
  #gráfico
  require(ggplot2)
  dados<-data.frame(x,yr,zr) #transformando os dados em data frame
  ggplot(dados,aes(x=x,y=yr))+
    geom_line(aes(col="Amostra fortificada"))+
    geom_line(aes(y=zr,col="amostra em branco"))
}
