pdf("saito_transform.pdf")
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
plot(foo$fdata.est[which(datanew$class=="Cyl"),], col="yellow",
     main = "Cylinder", xlab = "", ylab = "")
lines(foo$fdata.est[278,], lty=1,lwd=2)

plot(foo$fdata.est[which(datanew$class=="Bel"),], col="red",
     main = "Bell", xlab = "", ylab = "")
lines(foo$fdata.est[117 ,], lty=1,lwd=2)

plot(foo$fdata.est[which(datanew$class=="Fun"),], col="blue",
     main = "Funnel", xlab = "", ylab = "")
lines(foo$fdata.est[470,], lty=1,lwd=2)
dev.off()
