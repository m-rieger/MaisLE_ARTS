#### data preparation to loop through models ####
#################################################-

## 1) m1 to get detR for direct.ab
m1 <- data.frame(model = "m1", site = c("maisC", "maisD"), meth = "direct.ab", 
                 pos = "raw", AcSc = "ac1-sc1", detR = "all")
# ## 2) m2 to get AcSc for all methods
# m21 <- data.frame(model = "m2", site = c("maisC", "maisD"), meth = "direct.ab", 
#                   pos = "raw", AcSc = "ac1-sc1", detR = "CHOOSE")
# m22 <- data.frame(model = "m2", site = c("maisC", "maisD"), meth = "direct.in", 
#                   pos = "raw", AcSc = "ac4-sc2", detR = NA)
# m23 <- data.frame(model = "m2", site = "maisC", meth = "omni.ab", 
#                   pos = "raw", AcSc = "ac1-sc1", detR = NA)
# m24 <- data.frame(model = "m2", site = "maisC", meth = "omni.ml", 
#                   pos = "raw", AcSc = "ac2-sc2", detR = NA)
# ## 3) m3 to get mean or raw positions for all methods
# m31 <- data.frame(model = "m3", site = c("maisC", "maisD"), meth = "direct.ab", 
#                   pos = "all", AcSc = "CHOOSE", detR = "CHOOSE")
# m32 <- data.frame(model = "m3", site = c("maisC", "maisD"), meth = "direct.in", 
#                   pos = "all", AcSc = "CHOOSE", detR = NA)
# m33 <- data.frame(model = "m3", site = "maisC", meth = "omni.ab", 
#                   pos = "all", AcSc = "CHOOSE", detR = NA)
# m34 <- data.frame(model = "m3", site = "maisC", meth = "omni.ml", 
#                   pos = "all", AcSc = "CHOOSE", detR = NA)
## 4) m4 to get best method
m41 <- data.frame(model = "m4", site = c("maisC", "maisD"), meth = "direct.ab", 
                  pos = "CHOOSE", AcSc = "CHOOSE", detR = "CHOOSE")
m42 <- data.frame(model = "m4", site = c("maisC", "maisD"), meth = "direct.in", 
                  pos = "CHOOSE", AcSc = "CHOOSE", detR = NA)
m43 <- data.frame(model = "m4", site = "maisC", meth = "omni.ab", 
                  pos = "CHOOSE", AcSc = "CHOOSE", detR = NA)
m44 <- data.frame(model = "m4", site = "maisC", meth = "omni.ml", 
                  pos = "CHOOSE", AcSc = "CHOOSE", detR = NA)
## 5) m5 to model position accuracy
m5 <- data.frame(model = "m5", site = c("maisC", "maisD"), meth = "CHOOSE", 
                  pos = "CHOOSE", AcSc = "CHOOSE", detR = "CHOOSE")

## merge needed variables for model datasets
df.mod <- do.call(rbind, list(m1, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44, m5))

write.csv(df.mod, "./data/Lookup_model.csv", row.names = F)



