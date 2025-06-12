modeR <- rep(NA, 500)

for (i in 1:500) {
  modeR[i] <- Mode(unlist(lapply(R, '[[', i))[1:10000])
  if (is.na(modeR[i])) {
    modeR[i] <- mean(unlist(lapply(R, '[[', i))[1:10000])
  }
}
library(caret)
confusionMatrix(as.factor(modeR), as.factor(R_long_true))

acc <- rep(NA, 10000)
for (i in 1:10000) {
  acc[i] <- sum(R[[i]]==R_long_true)/500
}

Y_pred <- rep(NA, 500)

for (i in 1:500) {
  Ri <- Mode(unlist(lapply(R, '[[', i))[1:10000])
  if (is.na(Ri)) {
    Ri <- mean(unlist(lapply(R, '[[', i))[1:10000])
  }
  mui <- mean(unlist(lapply(mu, '[[', Ri))[1:10000])
  thetai <- mean(unlist(lapply(theta, function(x) x[i, Ri]))[1:10000])
  
  Y_pred[i] <- mui + thetai
}