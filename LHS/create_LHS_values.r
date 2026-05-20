# LHS ensemble

filename = paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.753/yelmox/LHS/lhs_np9_ns600.txt",sep="")
txt=read.table(filename,header=F, sep="")
a=as.matrix(txt)
colnames(a)=c("itmb","itmc","q","z0","cfngs","nff","enh","k","fp")    # order of columns depending on Yelmo output
nr=nrow(a)
nc=ncol(a)

# convert 0 to 1 numbers to variable domains
# itmc, cgrz, kppgrz, fp, bt0, cffrzn, cfstrm, enhshr

# itmb -2.4,-3.0
a[,1]=a[,1]*(-3.000-(-2.400))+(-2.400)

# itmc -24,-30
a[,2]=a[,2]*(-30.000-(-24.000))+(-24.000)

# q (0,1)
a[,3]=a[,3]*(1.000-(0.000))+(0.000)

# z0 (-600,0)
a[,4]=a[,4]*(0.000-(-600.00))+(-600.000)

# cfngs (0.1,0.5)
a[,5]=a[,5]*(0.500-(0.100))+(0.100)

# nffdlt (0.02,0.05)
a[,6]=a[,6]*(0.050-0.020)+0.020

#enh (1,4)
a[,7]=a[,7]*(4.000-1.000)+1.000

# k (1,10)
a[,8]=a[,8]*(10.000-1.000)+1.000

# fp (0.02, 0.1)
a[,9]=a[,9]*(0.100-0.020)+0.020

b=a
b[,1]=sprintf('%.3f', a[,1])
b[,2]=sprintf('%.3f', a[,2])
b[,3]=sprintf('%.3f', a[,3])
b[,4]=sprintf('%.3f', a[,4])
b[,5]=sprintf('%.3f', a[,5])
b[,6]=sprintf('%.3f', a[,6])
b[,7]=sprintf('%.3f', a[,7])
b[,8]=sprintf('%.3f', a[,8])
b[,9]=sprintf('%.3f', a[,9])

pairs(a, upper.panel=NULL)

# write the values in a txt file. 
out=format(b,digits=3,width=10,justify="right", scientific=F)
filename = paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.753/yelmox/LHS/lhs_np9_ns600_values.txt", sep="")
write.table(out,file=filename,row.names=FALSE,col.names=TRUE,quote=FALSE)

