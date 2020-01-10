cl <- apply(et_cl, 1, mean)
co <- apply(et_co, 1, mean)
xt <- xtable(cbind(cl, co), digits = 2)
rownames(xt) <- c('fun_1', 'fun_2', 'fun_3', 'ft1_1', 'ft1_2', 'ft1_3', 'ft2_1', 'ft2_2', 'ft2_3', 'ft3_1', 'ft3_2', 'ft3_3', 'ft4_1', 'ft4_2', 'ft4_3')

df <- cbind(cl, co)
df_split <- split(df, 1:3)
x1 <- df_split$'1'[1:5]
x2 <- df_split$'1'[6:10]
x3 <- df_split$'2'[1:5]
x4 <- df_split$'2'[6:10]
x5 <- df_split$'3'[1:5]
x6 <- df_split$'3'[6:10]

pdf()
par(mfrow=c(3,1))
plot(x1, type = 'b', main = 'Small dataset (n=20)', ylab = 'Computation time (sec)', xaxt='n', col = 'blue', lty = 1, pch = 19, xlab = 'Functions')
lines(x2, type = 'b', col = 'red', lty = 1, pch = 19)
axis(1, at = 1:5, labels = c('fun', 'ft1', 'ft2', 'ft3', 'ft4'))
legend("topright", legend = c("cluster", "coeff"), col = c('blue', 'red'), pch = 19, bty = "n", text.col = "black", horiz = F)

plot(x3, type = 'b', main = 'Medium dataset (n=100)', ylab = 'Computation time (sec)', xaxt='n', col = 'blue', lty = 1, pch = 19, xlab = 'Functions', ylim=c(0,100))
lines(x4, type = 'b', col = 'red', lty = 1, pch = 19)
axis(1, at = 1:5, labels = c('fun', 'ft1', 'ft2', 'ft3', 'ft4'))
legend("topright", legend = c("cluster", "coeff"), col = c('blue', 'red'), pch = 19, bty = "n", text.col = "black", horiz = F)

plot(x5, type = 'b', main = 'Large dataset (n=200)', ylab = 'Computation time (sec)', xaxt='n', col = 'blue', lty = 1, pch = 19, xlab = 'Functions', ylim = c(0,400))
lines(x6, type = 'b', col = 'red', lty = 1, pch = 19)
axis(1, at = 1:5, labels = c('fun', 'ft1', 'ft2', 'ft3', 'ft4'))
legend("topright", legend = c("cluster", "coeff"), col = c('blue', 'red'), pch = 19, bty = "n", text.col = "black")
dev.off()
