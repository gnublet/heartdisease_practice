a = c(5,2,4,6,1,3)
for (j in 2:length(a)) {
  key = a[j]
  i = j-1
  while (i>0 && a[i]>key){
    a[i+1] = a[i]
    i = i-1
  }
  a[i+1]= key
}

a
rm(a,i,j,key)
2:length(a)
