#'@export
#'@description
#'O presente pacote faz o tratamento de dados provenientes de uma análise
#cromatográfica informando quando a regressão linear, análise de premissas e
#valores extremos
#' @param x são as concentrações
#' @param y são as respostas cromatográficas
#'@example
#'x<- c(5.3,5.3,5.3,13.3,13.3,13.3,21.3,21.3,21.3,29.3,29.3,29.3,37.3,37.3,37.3,45.3,45.3,45.3)
#y<- c(423832,267896,370213,1545854,1093452,780549,2378730,2204902,2709421,3492212,3591964,3523445,4762476,4830668,4853425,5334327,5753813,5648307)
#Lin_UFVJM(x,y)

#Análise de premissas e linearidade
Lin_UFVJM<-function(x,y){
  library(car)
  #modelo de regressao
  md=lm(y~x)

  ANA<-summary(aov(y~factor(x)))
  # Montando a anova
  n<-length(y)
  nt<-length(levels(factor(x)))
  r<-n/nt
  trat<-factor(x)
  SQTotal<-sum(y^2) -(sum(y)^2)/n
  SQTotal
  b0<-coefficients(md)[1];b1<-coefficients(md)[2]
  SQReg<-b0*sum(y)+b1*sum(x*y) -(sum(y)^2)/n
  SQReg
  SQTrat<-sum(tapply(y,x, sum)^2)/r-(sum(y)^2)/n
  SQTrat
  SQRes<-SQTotal-SQTrat
  SQRes
  SQFa<-SQTrat-SQReg
  SQFa

  ####graus de liberdade
  glreg<-1
  gltotal<-n-1
  gltrat<-nt-1
  glres<-gltotal-gltrat
  glFa<-gltrat-glreg

  #### quadrados medios
  QMReg<-SQReg/glreg
  QMFa<-SQFa/glFa
  QMtrat<-SQTrat/gltrat
  QMres<-SQRes/glres
  GL<-c(glreg,glFa,gltrat,glres,gltotal)
  SQ<-round(c(SQReg,SQFa,SQTrat,SQRes,SQTotal),3)
  QM<-c(round(c(QMReg,QMFa,QMtrat,QMres),3),'')
  anava<-data.frame("GL"=GL,
                    "SQ"=SQ,
                    "QM"=QM,
                    "Fc"=c(round((QMReg/QMres),3),round((QMFa/QMres),3),
                           round((QMtrat/QMres),3),'',''),
                    "Pr>Fc"=c(round(1-pf(QMReg/QMres,glreg,glres),5),
                              round(1-pf(QMFa/QMres,glFa,glres),5),
                              round(1-pf(QMtrat/QMres,gltrat,glres),5),'',''))
  colnames(anava)[5]="Pr>Fc"
  rownames(anava)=c('Regressão','F.Ajuste','Tratamento','Resíduo','Total')
  ################################### ANOVA #######################

  cat('------------------------------------------------------------------------
Análise de variância do modelo linear\n------------------------------------------------------------------------\n')

  print(anava)

  ################################### Estatística t do modelo #######################
  tab2<-summary(md)
  TabReg<-tab2$coefficients
  row.names(TabReg)<-c('b0','b1')
  cat('------------------------------------------------------------------------
Estatística t do modelo\n------------------------------------------------------------------------\n')

  print(TabReg)
  cat("R² : ", tab2$r.squared, "\n")

  #################################################################################
  #Teste de normalidade
  pvalor.shapiro <- shapiro.test(md$residuals)$p.value
  cat("\n------------------------------------------------------------------------\nTeste de normalidade dos resíduos (Shapiro-Wilk)\n")
  cat("p-valor: ", pvalor.shapiro, "\n")
  if (pvalor.shapiro <= 0.05) {
    cat("ATENÇÃO: a 95% de confiança, os resíduos não podem ser considerados normais.\n------------------------------------------------------------------------\n")
  }
  if (pvalor.shapiro > 0.05) {
    cat("De acordo com o teste de Shapiro-Wilk a 95% de confiança, os resíduos podem ser considerados normais.\n------------------------------------------------------------------------\n\n")
  }

  #################################################################################
  #Teste da Homogenidade da variância Brown-Forsythe
  bf<-factor(x)
  t<-length(levels(bf))
  r<-as.numeric(table(bf))
  bf<-factor(x)
  zdados1<-matrix(0,length(y),1)
  rp<-0
  for(k in 1:length(y)) {
    zdados1[k]<-abs(y[k]-median(y[(rp+1):(rp+r[bf[k]])]))
    if(k<length(y)){if(x[k]!=x[k+1]){rp<-sum(r[1:bf[k]])}}
  }
  pvalor.bf<-summary(aov(zdados1 ~ bf))[[1]][1,5]

  cat('\n------------------------------------------------------------------------\nTeste de homogeneidade de variâncias (Brown-Forsythe) \n')
  cat('valor-p: ',pvalor.bf, '\n')
  if(pvalor.bf>0.05){cat('De acordo com o teste de Brown-Forsythe, a 95% de Confianca, as variâncias podem ser consideradas homogêneas!\n')
    cat('------------------------------------------------------------------------\n')}
  else{cat('ATENÇÃO: a 95% de confiança, as variâncias não podem ser consideradas homogêneas.\n')
    cat('------------------------------------------------------------------------\n')
  }

  #################################################################################
  #Teste de Independencia dos residuos Durbin-Watson
  pvalor.DW <-durbinWatsonTest(md)$p
  cat("\n------------------------------------------------------------------------\nTeste para independência dos resíduos(Durbin-Watson)\n")
  cat("p-valor: ", pvalor.DW, "\n")
  if (pvalor.DW <= 0.05) {
    cat("ATENÇÃo: a 95% de confiança,os resíduos são dependentes.\n------------------------------------------------------------------------\n")
  }
  if (pvalor.DW > 0.05) {
    cat("De acordo com o teste de Durbin-Watson a 95% de confiança, os resíduos são independentes.\n------------------------------------------------------------------------\n\n")
  }

  ############## Parte grafica
  par(mfrow=c(2,2))
  residuo<-md$residuals
  #valor predito
  yp <- predict.lm(md)
  #residuo Jaknife
  Jack <- rstudent(md)

  qq<-data.frame(y,Jack,yp,ni=1:n)
  dd<-qq[order(abs(qq[,2]), decreasing=TRUE),]
  pts<-ifelse(abs(dd[,2])>qt(0.975,n-2),1,0)*dd$ni
  pontos<-pts[min(which(pts!=0)):max(which(pts!=0))]
  soma=sum(pontos)

  #Grafico1
  #qq-plot
  qe<-function(r){
    n<-length(r)
    res<-c()
    for(i in 1:n){
      res[i]<-(i-3/8)/(n+1/4)
    }
    res
  }
  resord<-sort(residuo)
  resnorm<-qnorm(qe(residuo))
  #plot(resord,resnorm)

  ### Ajustado x preditos
  plot(md,which=1,ann=FALSE,sub="")
  title(xlab=" ", ylab=" ")
  ### QQplot
  plot(md,which=2,ann=FALSE,sub="")
  title(xlab=" ", ylab=" ")

  #Grafico 2
  plot(yp,Jack,ylim=c(-3,3),xlab='Valores Preditos',ylab='Resíduo de Pearson',main='',
  )
  abline(h=0)
  abline(h=qt(0.975,n-2),lty = 2)
  abline(h=qt(0.025,n-2), lty = 2)

  ### Distancia de Cook
  plot(md,which=4,ann=FALSE,sub="")
  title(xlab=" ", ylab=" ")

  #Pontos discrepantes via JK
  cat("Pontos discrepantes pela análise do resíduo Jackknife \n")
  cat('\n')
  if (soma == 0) {
    cat("Não há pontos discrepantes segundo o resído Jackknife.\n------------------------------------------------------------------------\n")
  }
  if (soma > 0) {
    cat("De acordo com o resíduo Jackknife, os pontos são:", pontos,"\n------------------------------------------------------------------------\n\n")
  }
  cat("Equação da reta:f(x)=",b1,"x +",b0," ")

  #Saida
  out<-list()
  out$anava
  out$pvalor.shapiro
  out$pvalor.levene
  #out$pontos
  invisible(out)

}
