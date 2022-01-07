library(Metrics)

Metrics::mape(c(1,2,3,4), c(2,2,2,2))

actual <- c(1, 1, 1, 0, 0, 0)
predicted <- c(0.9, 0.8, 0.4, 0.5, 0.3, 0.2)
Metrics::auc(actual, predicted)

df = fread('D:/Github/ACQ/ACQ.Excel/Tests/metrics.csv')

res = ldply(seq(nrow(df)), function(i){
  data.frame(auc1 = Metrics::auc(df$z[1:i], df$pred[1:i]), auc2=gbm.roc.area(df$z[1:i], df$pred[1:i]))
})

r <- rank(df$pred)
n_pos <- as.numeric(sum(df$z == 1))
n_neg <- length(df$z ) - n_pos
(sum(r[df$z == 1]) - n_pos * (n_pos + 1)/2)/(n_pos * n_neg)
sprintf('%.20f', Metrics::auc(df$z, df$pred))
sprintf('%.20f', (sum(r[df$z == 1])-1 - n_pos * (n_pos + 1)/2)/(n_pos * n_neg))
sprintf('%.20f', gbm.roc.area(df$z, df$pred))
sprintf('%.20f', binmodel_roc(df$z, df$pred))

plot_binmodel_cdf(df$z, df$pred)
binmodel_ks(df$z, df$pred)

mean(rank(df$pred)[df$z == 1]) * n_pos

sum(r[actual == 1]-1) 
169642