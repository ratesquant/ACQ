library(data.table)

df = fread('D:/Github/ACQ/ACQ.Excel/regression_data.csv')

summary( lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu, data = df) )
summary( lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu - 1, data = df) )
summary( lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu, data = df, weights = df$w) )

lm_wi = lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu-1, data = df, weights = df$w)

lm_w = lm(formula = y ~ age + sex + bmi + map + tc + ldl + hdl + tch + ltg + glu, data = df, weights = df$w)
yp = predict(lm_w, df)

cc(summary(lm_wi)$coefficients)

ssd <-function(x, w){
  xm = weighted.mean(x, w)
  return (sum( w * (x - xm)^2 ))
}

weighted.mean(yp, df$w)

mss = ssd(yp, df$w)
rss = ssd(df$y - yp, df$w)
mss / (mss + rss)
summary(lm_w)$r.squared