library(lhs)

#**************************************************************************************
# 3D LHS 
a=randomLHS(600,9)

a=round(a, digits=3)
#for(i in ncol(a)){
#  a[,i]=sprintf('%.3f', a[,i])
#}

out=format(a,digits=3,width=10,justify="right", scientific=F)
filename = paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.753/yelmox/LHS/lhs_np9_ns600.txt")
write.table(out,file=filename,row.names=FALSE,col.names=FALSE,quote=FALSE)
cat("LHS table written: ",filename,"\n")

#pairs(a,upper.panel=NULL, pch=20)

