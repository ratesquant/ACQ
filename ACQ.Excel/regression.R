library(data.table)

df = fread('regression_data.csv')

res = lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu, data = df)
res = lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu - 1, data = df)
res = lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu, data = df, weights = df$w)
res = lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu-1, data = df, weights = df$w)

summary(res)
AIC(res)
BIC(res)

yp = predict(res, df)

cc(summary(res)$coefficients)

ssd <-function(x, w){
  xm = weighted.mean(x, w)
  return (sum( w * (x - xm)^2 ))
}

ss <-function(x, w){
  return (sum( w * x^2 ))
}

weighted.mean(yp, df$w)


mss = ssd(yp, df$w)
rss = ssd(df$y - yp, df$w)
mss / (mss + rss)

mss = ss(yp, df$w)
rss = ss(df$y - yp, df$w)
mss / (mss + rss)
summary(res)$r.squared

