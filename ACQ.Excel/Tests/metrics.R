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
binmodel_ks_ex(df$z, df$pred)

binmodel_ks(c(0,1,1,0), c(0.2, 0.5, 0.6, 0.3))
binmodel_ks_ex(c(0,1,1,0), c(0.2, 0.5, 0.6, 0.3))

mean(rank(df$pred)[df$z == 1]) * n_pos

sum(r[actual == 1]-1) 
169642

x = c(0.50,	0.75,	1.00,	1.25,	1.50,	1.75,	1.75,	2.00,	2.25,	2.50,	2.75,	3.00, 	3.25, 	3.50, 	4.00, 	4.25, 	4.50, 	4.75, 	5.00, 	5.50)
y = c(0, 	0, 	0, 	0, 	0, 	0 ,	1, 	0, 	1, 	0, 	1, 	0, 	1, 	0 ,	1, 	1, 	1, 	1, 	1, 	1) 

model = glm(y ~ x, binomial(link = "logit"))
summary(model)
yp= predict(model, type = 'response')
binmodel_ks(y, yp)
plot_binmodel_predictions(y, yp)
gbm.roc.area(y, yp)

n1 <- sum(yp)
n <- length(yp)
((mean(rank(yp)[y > 0]) - (n1 + 1)/2)/(n - n1))


