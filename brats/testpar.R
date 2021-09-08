library('foreach')
library('doFuture')
registerDoFuture()

plan(sequential)
plan(multisession, workers = 6L)
a <- foreach(i=1:10)%dopar%{i}
print(a)