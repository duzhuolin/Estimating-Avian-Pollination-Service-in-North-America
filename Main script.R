#### jags function for SOCB model
library(data.table)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(jagsUI)
library(MCMCvis)
library(mgcv)
library(rjags)
library(stringr)
library(TruncatedNormal)
library(renv)

ind = "index"
lci = "lci"
uci = "uci"
year = "year"
base.yr = 1970
end.yr = 2017
popsource = "Pop.source"

fig_palette <- c("#D55E00", 
                 "#56B4E9", 
                 "#009E73", 
                 "#E69F00", 
                 "#0072B2", 
                 "#CC79A7", 
                 "#999999")

set.seed(2019)

pollinators = read.csv("input/avian nectarivore data for computing FVR.csv",
                       stringsAsFactors = F)
nectarivore = pollinators$species

popest = read.csv("input/Rosenberg et al species list.csv",
                  stringsAsFactors = F)
popest = popest[popest$species %in% nectarivore, ]
popest$popse = ((popest$popestuci-popest$popestlci)/(1.96*2))

indicesraw = read.csv("input/Rosenberg et al annual indices of abundance.csv",
                      stringsAsFactors = F)
indicesraw = indicesraw[order(indicesraw$species,indicesraw$year),]
indicesraw = indicesraw[indicesraw$species %in% nectarivore, ]


####################### GAM smoothing of population indices
sps = popest$species

for(ss in sps){
  wss = which(indicesraw$species == ss)
  tmp = indicesraw[wss,]
  
  fyr = unique(tmp$firstyear)
  lyr = unique(tmp$lastyear)
  
  torep = which(indicesraw$species == ss & (indicesraw$year >= fyr & indicesraw$year <= lyr))
  
  if(any(!is.na(tmp$lci.raw))){
    wprec = T
  }else{
    wprec = F
  }
  
  tmpd = tmp[which(!is.na(tmp$index.raw)),]
  
  
  
  nknots = min(c(11,max(floor(nrow(tmpd)/3),3)))
  if(ss %in% c("Cackling Goose","Greater White???fronted Goose")){nknots = 3} 
  
  
  
  form = as.formula(paste("lindex","~",
                          "s(year,k =",nknots,")"))
  
  
  
  
  ncounts = nrow(tmpd)
  ###### building the GAM basis function
  # gam basis functions created using function jagam from mgcv package          
  
  yminy = min(tmpd$year,na.rm = T)
  ymaxy = max(tmpd$year,na.rm = T)
  yearvec = tmpd$year-(yminy-1)
  yrs = seq(yminy,ymaxy,by = 1)
  nyears = length(yrs)
  ymin = min(yearvec)
  ymax = max(yearvec)
  tmpd$lindex = log(tmpd$index.raw)
  lindex = tmpd$lindex
  
  
  preddat = data.frame(lindex = 1,
                       year = yrs)
  
  
  
  if(wprec){
    ### if true then jags model runs 
    #### setting of the number of knots
    
    
    
    preci = 1/(((log(tmpd[,"uci.raw"])-log(tmpd[,"lci.raw"]))/(1.96*2))^2)
    
    
    
    gamprep = jagam(formula = form,
                    data = tmpd,
                    file = "models/tempgam.txt",
                    centred = T)
    
    gamprep.pred = jagam(formula = form,
                         data = preddat,
                         file = "models/tempgampred.txt",
                         centred = T)
    
    
    
    dat = list(X = gamprep$jags.data$X,
               S1 = gamprep$jags.data$S1,
               ncounts = nrow(tmpd),
               lindex = lindex,
               nknots = nknots,
               preci = preci,
               X.pred = gamprep.pred$jags.data$X,
               zero = gamprep$jags.data$zero)
    
    mg <- jags.model(data = dat,
                     file = paste0("models/GAM model smoothing indices jagam.txt"),
                     n.chains = 3)
    
    adaptest <- adapt(object = mg,
                      n.iter = 10)
    
    while(adaptest == F){
      adaptest <- adapt(object = mg,
                        n.iter = 1000)
      
    }
    
    nburn = 10000
    mgo = coda.samples(mg,
                       c("ind.pred"
                         #"rho"
                         #"mu",
                         #"b"
                       ),
                       n.iter = 10000,
                       n.burnin = nburn,
                       thin = 10)
    
    
    mgosum = summary(mgo)
    
    predg = mgosum$quantiles
    
    predgs = mgosum$statistics
    
    
    indicesraw[torep,"index"] <- exp(predg[grep(row.names(predg),pattern = "ind.pred"),"50%"])
    indicesraw[torep,"lci"] <- exp(predg[grep(row.names(predg),pattern = "ind.pred"),"2.5%"])
    indicesraw[torep,"uci"] <- exp(predg[grep(row.names(predg),pattern = "ind.pred"),"97.5%"])
    
    
    
  }else{ ### else if wprec
    
    
    m1 = gam(formula = form,
             data = tmpd)
    
    
    pred = predict(m1,newdata = preddat,
                   type = "link",
                   se.fit = T)
    
    
    
    
    indicesraw[torep,"index"] <- exp(pred$fit)
    indicesraw[torep,"lci"] <- exp(pred$fit-(1.96*pred$se.fit))
    indicesraw[torep,"uci"] <- exp(pred$fit+(1.96*pred$se.fit))
    
    
  } ### end if wprec
  
  
}

#save(indicesraw,file = "temp_output/post GAM indices.RDATA")
write.csv(indicesraw, "temp_output/post_GAM_indices_nectarivore.csv", row.names = FALSE)

######################### end GAM smoothing of annual indices


indices = indicesraw


yrs = sort(unique(indices$year))


indices[,"se"] <- ((indices[,uci]-indices[,lci])/(1.96*2))

indices <- indices[order(indices$firstyear,indices$species,indices$year),]

splist <- unique(indices[,c("species","firstyear","lastyear")])
splist3 = merge(splist,popest,by.x = "species",by.y = "species")

#adding taxa (order), habitat, and specialization level to splist3
match_idx <- match(splist3$species, pollinators$species)
splist3$specialization_lvl <- pollinators$specialization_lvl[match_idx]
splist3$Order <- pollinators$Order[match_idx]
splist3$Habitat <- pollinators$Habitat[match_idx]

splist3$spfactor <- factor(splist3$species,
                           levels = splist3$species,
                           ordered = T)

splist3$spfact <- as.integer(splist3$spfactor) 

indices$spfactor = factor(indices$species,
                          levels = splist3$species,
                          ordered = T)


splist = splist3
base.i <- rep(NA, length = length(unique(indices$species)))
names(base.i) <- unique(indices$species)
base.se.sp <- rep(NA, length = length(unique(indices$species)))
names(base.se.sp) <- unique(indices$species)
se = "se"

for (sn in 1:nrow(splist)) {
  
  s = as.character(splist[sn,"species"])
  
  
  
  r <- which(indices[,"species"] == s)
  rmpre = which(indices[,"species"] == s &
                  is.na(indices[,ind]) &
                  indices[,"year"] < splist[sn,"firstyear"])
  rmpost = which(indices[,"species"] == s &
                   is.na(indices[,"se"]) &
                   indices[,"year"] > splist[sn,"lastyear"])
  
  
  base.s <- indices[which(indices$species == s & indices$year == base.yr),ind] #stores the base index value 
  base.se <- indices[which(indices$species == s & indices$year == base.yr),se] #stores the base se value 
  if(is.na(base.s)){
    byr <- min(indices[which(indices$species == s & !is.na(indices[,ind])),"year"],na.rm = T)
    base.s <- indices[which(indices$species == s & indices$year == byr),ind] #stores the base index value 
    base.se <- indices[which(indices$species == s & indices$year == byr),se] #stores the base se value 
    
  }
  for (y in r) {
    
    indices[y,"index.s"] <- (indices[y,ind])/base.s # standardized index WRT base year
    indices[y,"cvar.s"] <- (((indices[y,se]^2)/(indices[y,ind]^2))+((base.se^2)/(base.s^2))) 
    indices[y,"logthetahat"] <- log(indices[y,"index.s"])      
    indices[y,"prec.logthetahat"] <- 1/log(1+(indices[y,"cvar.s"]))
  }
  print(s)
}

indices = merge(indices,splist[,c("species","spfactor","spfact")],by = "species")
indices$year_i = (indices$year - base.yr)+1

indices = indices[order(indices$spfact,indices$year_i),]




splist = splist[order(splist$spfact),]
nspecies <- nrow(splist)
nyears = max(indices$year_i)

splist$g1 = as.integer(factor(splist$Winter.Biome))
splist$g2 = as.integer(factor(splist$Breeding.Biome))

grps = as.matrix(splist[,c("g1","g2")])

ngroups1 = max(grps[,1])
ngroups2 = max(grps[,2])

ident = function(x){
  return(x)
}
logthetahat <- as.matrix(tapply(indices[,"logthetahat"],indices[,c("spfact","year_i")],ident))
prec.logthetahat <- as.matrix(tapply(indices[,"prec.logthetahat"],indices[,c("spfact","year_i")],ident))

ne = splist$popest #population estimate
tau.ne = 1/splist$popse^2 #precision of population estimate



# 
# 
yest = (splist$year_est1-base.yr)+1 # first year over which species population should be averaged

yavg = 1+(splist$year_est2-splist$year_est1) ## number of years over which to average the species population estimate





##### species group indexing (breeding-wintering biome combination)

subgrps = unique(grps)
subgrps = subgrps[order(subgrps[,1],subgrps[,2]),]
nsubgrps = table(subgrps[,1])
subgrpsl = matrix(NA,ncol = ngroups1,nrow = max(nsubgrps))

nsppsubbiomesmat = table(grps[,1],grps[,2])

spsubgrpmat = array(NA,dim = c(ngroups1,ngroups2,max(nsppsubbiomesmat)))


for(i in 1:ngroups1){
  subs = subgrps[which(subgrps[,1] == i),2]
  
  subgrpsl[1:nsubgrps[i],i] = subs 
  for(j in subs){
    spsubgrpmat[i,j,1:nsppsubbiomesmat[i,j]] = which(grps[,1] == i & grps[,2] == j)
    
  }
}#i


#### breeding biome indexing

nsppbiomes = table(grps[,2])
spinbiomes = matrix(NA,nrow = max(nsppbiomes),ncol = max(ngroups2))
for(g in 1:ngroups2){
  
  spinbiomes[1:nsppbiomes[g],g] <- which(grps[,2] == g)
}


### wintering biome indexing
nsppwinters = table(grps[,1])
spinwinters = matrix(NA,nrow = max(nsppwinters),ncol = max(ngroups1))
for(g in 1:ngroups1){
  
  spinwinters[1:nsppwinters[g],g] <- which(grps[,1] == g)
}

### specialization level indexing

splist$spclfact = (factor(splist$specialization_lvl))
splist$spclfactn = as.integer(factor(splist$specialization_lvl))

nsppspcls = table(splist$spclfact)
spinspcls = matrix(NA,nrow = max(nsppspcls),ncol = length(nsppspcls))
for(s in 1:length(nsppspcls)){
  sn = names(nsppspcls)[s]
  spinspcls[1:nsppspcls[s],s] <- which(splist$specialization_lvl == sn)
}
nspcls = length(nsppspcls)
spcls = unique(splist[,c("specialization_lvl","spclfact","spclfactn")])
spcls = spcls[order(spcls$spclfactn),]

### habitat indexing

splist$hbfact = (factor(splist$Habitat))
splist$hbfactn = as.integer(factor(splist$Habitat))

nspphbs = table(splist$hbfact)
spinhbs = matrix(NA,nrow = max(nspphbs),ncol = length(nspphbs))
for(h in 1:length(nspphbs)){
  hn = names(nspphbs)[h]
  spinhbs[1:nspphbs[h],h] <- which(splist$Habitat == hn)
}
nhbs = length(nspphbs)
hbs = unique(splist[,c("Habitat","hbfact","hbfactn")])
hbs = hbs[order(hbs$hbfactn),]

### indexing for family summaries

splist$famfact = (factor(splist$Family))
splist$famfactn = as.integer(factor(splist$Family))

nsppfams = table(splist$famfact)
spinfams = matrix(NA,nrow = max(nsppfams),ncol = length(nsppfams))
for(f in 1:length(nsppfams)){
  fn = names(nsppfams)[f]
  spinfams[1:nsppfams[f],f] <- which(splist$Family == fn)
}
nfams = length(nsppfams)
fams = unique(splist[,c("Family","famfact","famfactn")])
fams = fams[order(fams$famfactn),]

### indexing for order summaries

splist$ordfact = (factor(splist$Order))
splist$ordfactn = as.integer(factor(splist$Order))

nsppords = table(splist$ordfact)
spinords = matrix(NA,nrow = max(nsppords),ncol = length(nsppords))
for(o in 1:length(nsppords)){
  on = names(nsppords)[o]
  spinords[1:nsppords[o],o] <- which(splist$Order == on)
}
nords = length(nsppords)
ords = unique(splist[,c("Order","ordfact","ordfactn")])
ords = ords[order(ords$ordfactn),]

#indexing for bird.group summaries

splist$birdgroupfact = (factor(splist$bird.group))
splist$birdgroupfactn = as.integer(factor(splist$bird.group))

nsppbirdgroups = table(splist$birdgroupfact)
spinbirdgroups = matrix(NA,nrow = max(nsppbirdgroups),ncol = length(nsppbirdgroups))
for(f in 1:length(nsppbirdgroups)){
  fn = names(nsppbirdgroups)[f]
  spinbirdgroups[1:nsppbirdgroups[f],f] <- which(splist$bird.group == fn)
}
nbirdgroups = length(nsppbirdgroups)
birdgroups = unique(splist[,c("bird.group","birdgroupfact","birdgroupfactn")])
birdgroups = birdgroups[order(birdgroups$birdgroupfactn),]


#indexing for migration summaries

splist$migratefact = (factor(splist$Migrate))
splist$migratefactn = as.integer(factor(splist$Migrate))

nsppmigrates = table(splist$migratefact)
spinmigrates = matrix(NA,nrow = max(nsppmigrates),ncol = length(nsppmigrates))
for(f in 1:length(nsppmigrates)){
  fn = names(nsppmigrates)[f]
  spinmigrates[1:nsppmigrates[f],f] <- which(splist$Migrate == fn)
}
nmigrates = length(nsppmigrates)
migrates = unique(splist[,c("Migrate","migratefact","migratefactn")])
migrates = migrates[order(migrates$migratefactn),]

#indexing for ai summaries

splist$aifact = (factor(splist$AI))
splist$aifactn = as.integer(factor(splist$AI))

nsppais = table(splist$aifact)
spinais = matrix(NA,nrow = max(nsppais),ncol = length(nsppais))
for(f in 1:length(nsppais)){
  fn = names(nsppais)[f]
  spinais[1:nsppais[f],f] <- which(splist$AI == fn)
}
nais = length(nsppais)
ais = unique(splist[,c("AI","aifact","aifactn")])
ais = ais[order(ais$aifactn),]

#indexing for native summaries

splist$nativefact = (factor(splist$native))
splist$nativefactn = as.integer(factor(splist$native))

nsppnatives = table(splist$nativefact)
spinnatives = matrix(NA,nrow = max(nsppnatives),ncol = length(nsppnatives))
for(f in 1:length(nsppnatives)){
  fn = names(nsppnatives)[f]
  spinnatives[1:nsppnatives[f],f] <- which(splist$native == fn)
}
nnatives = length(nsppnatives)
natives = unique(splist[,c("native","nativefact","nativefactn")])
natives = natives[order(natives$nativefactn),]



#imputing the missing data with assumptions of no-change from most recent year with data and gradually decreasing precision

wspecieslate = as.integer(splist[which(splist$firstyear > 1970),"spfact"])
yearswo1late = rep(1,nspecies)
yearswo2late = (splist[,"firstyear"])-1970


wspeciesearly = as.integer(splist[which(splist$lastyear < 2017),"spfact"])
yearswo1early = (splist[,"lastyear"]+1)-1969
yearswo2early = rep(nyears,nspecies)

yearsw1 = (splist[,"firstyear"])-1969
yearsw2 = (splist[,"lastyear"])-1969



prec.powdrop = 2 #exponential function for decreasing precision with years (i.e., precision decreases with square of the years since real data)

for( s in wspecieslate) { # wspecieslate = vector of species that don't have data in year-1
  for(y in yearswo1late[s]:yearswo2late[s]){
    
    logthetahat[s,y] <- logthetahat[s,yearsw1[s]]
    prec.logthetahat[s,y] <- prec.logthetahat[s,yearsw1[s]]/((yearsw1[s]-y)^prec.powdrop)# 
  }}


### imputing missing data for species without data at the end of the time series
## currently assumes that the precision decreases with the square of the number of years since the last data
for( s in wspeciesearly) { # wspecieslate = vector of species that don't have data in year-1
  for(y in yearswo1early[s]:yearswo2early[s]){
    
    logthetahat[s,y] <- logthetahat[s,yearsw2[s]]
    prec.logthetahat[s,y] <- prec.logthetahat[s,yearsw2[s]]/((y-yearsw2[s])^prec.powdrop)# 
    
    
    
  }}

############################# MCMC sampling

data.jags = list(nspecies = nspecies,
                 nyears = nyears,
                 logthetahat = logthetahat,
                 prec.logthetahat = prec.logthetahat,
                 ne = ne,
                 tau.ne = tau.ne,
                 yest = yest,
                 ngroups1 = ngroups1,
                 ngroups2 = ngroups2,
                 grps = grps,
                 yavg = yavg,
                 
                 spinbiomes = spinbiomes,
                 nsppbiomes = nsppbiomes,
                 
                 spinspcls = spinspcls,
                 nsppspcls = nsppspcls,
                 nspcls = nspcls,
                 
                 spinhbs = spinhbs,
                 nspphbs = nspphbs,
                 nhbs = nhbs,
                 
                 spinfams = spinfams,
                 nsppfams = nsppfams,
                 nfams = nfams,
                 
                 spinords = spinords,
                 nsppords = nsppords,
                 nords = nords,
                 
                 spinbirdgroups = spinbirdgroups,
                 nsppbirdgroups = nsppbirdgroups,
                 nbirdgroups = nbirdgroups,
                 
                 spinmigrates = spinmigrates,
                 nsppmigrates = nsppmigrates,
                 nmigrates = nmigrates,
                 
                 spinais = spinais,
                 nsppais = nsppais,
                 nais = nais,
                 
                 spinwinters = spinwinters,
                 nsppwinters = nsppwinters,
                 
                 spinnatives = spinnatives,
                 nsppnatives = nsppnatives,
                 nnatives = nnatives,
                 
                 subgrpsl = subgrpsl,
                 nsubgrps = nsubgrps,
                 nsppsubbiomesmat = nsppsubbiomesmat,
                 spsubgrpmat = spsubgrpmat)

params = c("expmu1",
           #"expm2",
           "expmu2",
           "Psi",
           "tau",
           "NEst",
           "Nsum",
           "N",
           "Nlost",
           "Nalost",
           "Nlost.S",
           "Nlost.fam",
           "Nlost.ord",
           "Nlost.spcl",
           "Nlost.hb",
           "Nlost.biome",
           "Nalost.fam",
           "Nalost.ord",
           "Nalost.spcl",
           "Nalost.hb",
           "Nlost.migrate",
           "Nalost.migrate",
           "plost.migrate",
           "Nlost.birdgroup",
           "Nalost.birdgroup",
           "plost.birdgroup",
           "Nlost.ai",
           "Nalost.ai",
           "plost.ai",
           "plost.native",
           "Nlost.native",
           "Nalost.native",
           "plost.winter",
           "Nlost.winter",
           "Nalost.winter",
           "Nalost.biome",
           "plost.fam",
           "plost.ord",
           "plost.biome",
           "plost.spcl",
           "plost.hb",
           "plost",
           "plost.S",
           "Nsum.subgrp")

mod = "models/Rosenberg et al model.txt"



adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 30000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=10000           # Total number of steps to save per chain.
thinSteps=50                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps ) + burnInSteps) / nChains ) # Steps per chain.

t1 = Sys.time()


jagsMod = jags(data = data.jags,
               model.file = mod,
               n.chains = nChains,
               n.adapt = adaptSteps,
               n.burnin = burnInSteps,
               n.thin = thinSteps,
               n.iter = nIter+burnInSteps,
               parameters.to.save = params,
               parallel = T)

############################################ end MCMC sampling

#save.image(file = "temp_output/Full_Project_Snapshot_after_MCMC.RData") #save the entire working environment

load("temp_output/Full_Project_Snapshot_after_MCMC.RData")

q90 = function(x){
  quantile(x,probs = c(0.025,0.25,0.75,0.975))
}
sumq = MCMCsummary(jagsMod$samples,func = q90,Rhat = F,n.eff = F,func_name = c("lci","lqrt","uqrt","uci"))


sumq = data.frame(sumq)
names(sumq) <- c("mean","sd","lci95","med","uci95","lci","lqrt","uqrt","uci")
write.csv(sumq,paste0("temp_output/pollinator population change parameters NA loss.csv"))
#write.csv(sumqalt,paste0("population change parameters NA loss w neff.csv"))


ordl = as.data.frame(sumq[paste0("Nlost.ord[",1:nords,"]"),])
ordp = as.data.frame(sumq[paste0("plost.ord[",1:nords,"]"),])
names(ordp) = paste0("plost.ord.",names(ordp))
ordl = cbind(ordl,ordp)
ordl = cbind(ordl,ords)

spcll = as.data.frame(sumq[paste0("Nlost.spcl[",1:nspcls,"]"),])
spclp = as.data.frame(sumq[paste0("plost.spcl[",1:nspcls,"]"),])
names(spclp) = paste0("plost.spcl.",names(spclp))
spcll = cbind(spcll,spclp)
spcll = cbind(spcll,spcls)

hbl = as.data.frame(sumq[paste0("Nlost.hb[",1:nhbs,"]"),])
hbp = as.data.frame(sumq[paste0("plost.hb[",1:nhbs,"]"),])
names(hbp) = paste0("plost.hb.",names(hbp))
hbl = cbind(hbl,hbp)
hbl = cbind(hbl,hbs)

biomes = unique(splist[,c("g2","Breeding.Biome")])
biomes = biomes[order(biomes$g2),]
biomel = as.data.frame(sumq[paste0("Nlost.biome[",1:ngroups2,"]"),])
biomep = as.data.frame(sumq[paste0("plost.biome[",1:ngroups2,"]"),])
names(biomep) = paste0("plost.biome.",names(biomep))
biomel = cbind(biomel,biomep)
biomel = cbind(biomel,biomes)



alll = data.frame(sumq[c("Nlost","plost"),])
allp = data.frame(sumq[c("plost","Nlost"),])
names(allp) = paste0("plost",names(allp))
alll = cbind(alll,allp)
alll = alll[1,]
alll[,c(19:21)] <- NA
biomel[,21] <- NA
names(alll)[20] = "Group"
alll$nspecies = nspecies
biomel$nspecies = nsppbiomes
ordl$nspecies = nsppords
spcll$nspecies = nsppspcls
hbl$nspecies = nspphbs
nms = names(alll)

names(biomel) = nms
names(ordl) = nms
names(spcll) = nms
names(hbl) = nms

allsums = rbind(alll,
                biomel,
                ordl,
                spcll, 
                hbl, stringsAsFactors = FALSE)

allsout = allsums[,c("Group",
                     "nspecies",
                     "lci",
                     "med",
                     "uci",
                     "plostlci",
                     "plostmed",
                     "plostuci")]
for(j in c("lci","med","uci")){
  allsout[,j] <- signif(allsout[,j]/1e6,5)
}
for(j in c("plostlci","plostmed","plostuci")){
  allsout[,j] <- signif(allsout[,j],3)
}



#### biome summaries

biomepop = expand.grid(biome = 1:ngroups2,
                       yrs = 1:nyears)
biomepop$param = paste0("Nalost.biome[",biomepop$biome,",",biomepop$yrs,"]")

biomepopt = as.data.frame(sumq[biomepop$param,])
biomepop = cbind(biomepop,biomepopt)
biomepop = merge(biomepop,biomes,by.x = "biome",by.y = "g2")
biomepop$year = biomepop$yrs + (base.yr-1)
biomepop = biomepop[order(biomepop$Breeding.Biome,biomepop$year),]
biomlab = biomepop[which(biomepop$year == 2017),]

biomlab$label_x <- NA
biomlab$label_y <- NA

biomlab$label_x[biomlab$Breeding.Biome == "Eastern Forest"] <- 2015
biomlab$label_y[biomlab$Breeding.Biome == "Eastern Forest"] <- 7

biomlab$label_x[biomlab$Breeding.Biome == "Aridlands"] <- 2015
biomlab$label_y[biomlab$Breeding.Biome == "Aridlands"] <- -9

biomlab$label_x[biomlab$Breeding.Biome == "Western Forest"] <- 2015
biomlab$label_y[biomlab$Breeding.Biome == "Western Forest"] <- -43

pmain = ggplot(data = biomepop,aes(x = year,y = med/1e6))+
  geom_ribbon(aes(x = year,ymin = lci/1e6,ymax = uci/1e6,group = Breeding.Biome,fill = Breeding.Biome),alpha = 0.2)+
  geom_line(aes(colour = Breeding.Biome))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_label(
    data = biomlab, 
    aes(x = label_x, y = label_y, label = Breeding.Biome, colour = Breeding.Biome),
    fill = "white",
    label.size = NA, 
    hjust = 1,          
    size = 9,
    inherit.aes = FALSE 
  ) +
  labs(x = "Year",y = "Change in no. of avian nectarivores (Millions)")+
  scale_x_continuous(limits = c(1970, 2020), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = pretty, expand = c(0, 0.5)) +
  scale_fill_manual(values = fig_palette) + 
  scale_color_manual(values = fig_palette)+
  theme_minimal()+
  theme(legend.position = "none", panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 25, margin = margin(t = 5)),
        axis.text.y = element_text(color = "black", size = 25, margin = margin(r = 5)),
        axis.title = element_text(size = 25),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
  )

pdf(paste0("output/Breeding biome level avian pollinator population change trajectory loss.pdf"), width = 10, height = 7)
print(pmain)
dev.off()
png("output/Breeding biome level avian pollinator population change trajectory loss.png", width = 3000, height = 2500, res = 300)
print(pmain)
dev.off()



#### specialization level summaries

spclpop = expand.grid(spcl = 1:nords,
                      yrs = 1:nyears)
spclpop$param = paste0("Nalost.spcl[",spclpop$spcl,",",spclpop$yrs,"]")

spclpopt = as.data.frame(sumq[spclpop$param,])
spclpop = cbind(spclpop,spclpopt)
spclpop = merge(spclpop,spcls,by.x = "spcl",by.y = "spclfactn")
spclpop$year = spclpop$yrs + (base.yr-1)
spclpop = spclpop[order(spclpop$specialization_lvl,spclpop$year),]
spcllab = spclpop[which(spclpop$year == 2017),]

spcllab$label_x <- NA
spcllab$label_y <- NA

spcllab$label_x[spcllab$specialization_lvl == "specialized"] <- 2010
spcllab$label_y[spcllab$specialization_lvl == "specialized"] <- -12

spcllab$label_x[spcllab$specialization_lvl == "generalized"] <- 2010
spcllab$label_y[spcllab$specialization_lvl == "generalized"] <- -23

pmain = ggplot(data = spclpop,aes(x = year,y = med/1e6))+
  geom_ribbon(aes(x = year,ymin = lci/1e6,ymax = uci/1e6,group = specialization_lvl,fill = specialization_lvl),alpha = 0.2)+
  geom_line(aes(colour = specialization_lvl))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_label(data = spcllab, aes(x = label_x, y = label_y, label = specialization_lvl, colour = specialization_lvl),
             fill = "white",
             label.size = NA,  
             hjust = 1,        
             size = 9) +
  labs(x = "Year",y = "Change in no. of avian nectarivores (Millions)")+
  scale_x_continuous(limits = c(1970, 2020), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = pretty, expand = c(0, 0.5)) +
  scale_fill_manual(values = fig_palette) + 
  scale_color_manual(values = fig_palette)+
  theme_minimal()+
  theme(legend.position = "none", panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 25, margin = margin(t = 5)),
        axis.text.y = element_text(color = "black", size = 25, margin = margin(r = 5)),
        axis.title = element_text(size = 25),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
  )

pdf(paste0("output/specialization_lvl level avian pollinator population change trajectory loss.pdf"), width = 10, height = 7)
print(pmain)
dev.off()
png("output/specialization_lvl level avian pollinator population change trajectory loss.png", width = 3000, height = 2500, res = 300)
print(pmain)
dev.off()


#### habitat summaries

hbpop = expand.grid(hb = 1:nords,
                    yrs = 1:nyears)
hbpop$param = paste0("Nalost.hb[",hbpop$hb,",",hbpop$yrs,"]")

hbpopt = as.data.frame(sumq[hbpop$param,])
hbpop = cbind(hbpop,hbpopt)
hbpop = merge(hbpop,hbs,by.x = "hb",by.y = "hbfactn")
hbpop$year = hbpop$yrs + (base.yr-1)
hbpop = hbpop[order(hbpop$Habitat,hbpop$year),]
hblab = hbpop[which(hbpop$year == 2017),]

biomlab$label_x <- NA
biomlab$label_y <- NA

hblab$label_x[hblab$Habitat == "Forest"] <- 2010
hblab$label_y[hblab$Habitat == "Forest"] <- -53

hblab$label_x[hblab$Habitat == "Shrubland"] <- 2010
hblab$label_y[hblab$Habitat == "Shrubland"] <- -17

hblab$label_x[hblab$Habitat == "Woodland"] <- 2010
hblab$label_y[hblab$Habitat == "Woodland"] <- 25

pmain = ggplot(data = hbpop,aes(x = year,y = med/1e6))+
  geom_ribbon(aes(x = year,ymin = lci/1e6,ymax = uci/1e6,group = Habitat,fill = Habitat),alpha = 0.2)+
  geom_line(aes(colour = Habitat))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_label(data = hblab, aes(x = label_x, y = label_y, label = Habitat, colour = Habitat),
             fill = "white",
             label.size = NA,  
             hjust = 1,        
             size = 9) +
  labs(x = "Year",y = "Change in no. of avian nectarivores (Millions)")+
  scale_x_continuous(limits = c(1970, 2020), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = pretty, expand = c(0, 0.5)) +
  scale_fill_manual(values = fig_palette) + 
  scale_color_manual(values = fig_palette)+
  theme_minimal()+
  theme(legend.position = "none", panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 25, margin = margin(t = 5)),
        axis.text.y = element_text(color = "black", size = 25, margin = margin(r = 5)),
        axis.title = element_text(size = 25),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
  )

pdf(paste0("output/Habitat level avian pollinator population change trajectory loss.pdf"), width = 10, height = 7)
print(pmain)
dev.off()
png("output/Habitat level avian pollinator population change trajectory loss.png", width = 3000, height = 2500, res = 300)
print(pmain)
dev.off()


#### overall summaries

totp = as.data.frame(sumq[paste0("Nsum[",1:nyears,"]"),c("med","lci","uci")])
names(totp) = paste0("N_",names(totp))
totp$year_i = 1:nyears
totp$year = totp$year_i + (base.yr-1)


lossa = sumq[paste0("Nalost[",1:nyears,"]"),c("med","lci","uci")]
names(lossa) = paste0("Loss_",names(lossa))
lossa = cbind(totp,lossa)

lost = as.data.frame(sumq["Nlost",])

lost.s = as.data.frame(sumq[paste0("Nlost.S[",1:nspecies,"]"),c("med","lci","uci","lqrt","uqrt")])
names(lost.s) = paste0("Loss_",names(lost.s))
for(i in 1:ncol(lost.s)){
  lost.s[,i] <- signif(lost.s[,i],3)
}

lost.s = cbind(splist,lost.s)


# June 2021 - adding species level proportional loss ----------------------

plost.s = as.data.frame(sumq[paste0("plost.S[",1:nspecies,"]"),c("med","lci","uci","lqrt","uqrt")])
names(plost.s) = paste0("Proportional_loss_",names(plost.s))
for(i in 1:ncol(plost.s)){
  plost.s[,i] <- signif(plost.s[,i],3)
}
lost.s = cbind(lost.s,plost.s)


# June 2021 - adding starting and ending population sizes -----------------

n1.S = as.data.frame(sumq[paste0("N[",1:nspecies,",1]"),c("med","lci","uci","lqrt","uqrt")])
names(n1.S) = paste0("Estimated_population_size_1970_",names(n1.S))
for(i in 1:ncol(n1.S)){
  n1.S[,i] <- signif(n1.S[,i],3)
}
n2.S = as.data.frame(sumq[paste0("N[",1:nspecies,",",nyears,"]"),c("med","lci","uci","lqrt","uqrt")])
names(n2.S) = paste0("Estimated_population_size_2017_",names(n2.S))
for(i in 1:ncol(n2.S)){
  n2.S[,i] <- signif(n2.S[,i],3)
}

lost.s = cbind(lost.s,n1.S)
lost.s = cbind(lost.s,n2.S)

write.csv(lost.s,paste0("output/Avian pollinator cumulative population change by species.csv"))


# identify the most threatened species
lost_most <- as.data.table(lost.s)
high_threat <- lost_most[Proportional_loss_lci > 0 & Proportional_loss_med > 0.5]
fwrite(high_threat, "output/Avian pollinators with the greatest decline.csv", row.names = F)


#summarize th eproportion of declining species for tables 1 and S2
lost.s$decline = F
lost.s[which(lost.s$Loss_med > 0),"decline"] = T
allsout$node = row.names(allsout)

for(j in 1:nrow(allsout)){
  if(j == 1){
    S = lost.s[,"decline"]
    tf = table(S)/length(S)
    
  }else{
    nn = gsub(gsub(allsout[j,"node"],pattern = "Nlost.",fixed = T,replacement = ""),pattern = "\\[.*",replacement = "")
    nn = paste0(nn,"fact")
    
    if(grepl(nn,pattern = "biome")){
      nn = "Breeding.Biome" 
    }
    if(grepl(nn,pattern = "winter")){
      nn = "Winter.Biome" 
    }
    S = lost.s[which(lost.s[,nn] == allsout[j,"Group"]),"decline"]
    tf = table(S)/length(S)
  }
  allsout[j,"proportion.species.decline"] = round(as.numeric(tf["TRUE"]), 1)
  allsout[j,"n.species.decline"] = sum(S)
}
write.csv(allsout,"output/Summary-Net change in abundance across the North American avian pollinators.csv",  row.names = F)


#continental pollinator population trajectory
lostpy = 0.4*max(totp$N_med)
if(lost$med > 0){
  lostp = paste0(signif(lost$med/1e6,4)," Million avian pollinators lost [",signif(lost$lci/1e6,4),"-",signif(lost$uci/1e6,4),"]")
}else{
  lostp = paste0(signif(abs(lost$med/1e6),4)," Million avian pollinators gained [",signif(abs(lost$lci)/1e6,4),"-",signif(abs(lost$uci)/1e6,4),"]")
  
}
sppop = sumq[paste0("N[",rep(1:nspecies,each = nyears),",",rep(1:nyears,times = nspecies),"]"),]
sppop = as.data.frame(sppop)
sppop$spfact = rep(1:nspecies,each = nyears)
sppop$year_i = rep(1:nyears,times = nspecies)
sppop$year =  sppop$year_i + (base.yr-1)
sppop = merge(sppop,splist,
              by= "spfact")


spstartpop = sppop[which(sppop$year == base.yr),]
spstartpop$meanpopstart = spstartpop$mean
indices2 = merge(indices,spstartpop[,c("species","meanpopstart")],by = "species")

for(s in unique(indices2$species)){
  wsp = which(indices2$species == s)
  
  by = unique(indices2[wsp,"firstyear"])
  
  wspb = which(indices2$species == s & indices2$year == by)
  
  indices2[wsp,"index.s.raw"] <- indices2[wsp,"index.raw"]/indices2[wspb,"index.raw"]
  
}

indices2$rescindex = indices2$index.s*indices2$meanpopstart
indices2$rescindex.raw = indices2$index.s.raw*indices2$meanpopstart
sppop2 = merge(sppop,indices2[,c("species","rescindex","rescindex.raw","year")],
               by = c("species","year"))


for(j in c("N_med","N_lci","N_uci")){
  totp[,j] = totp[,j]/1e6
}

splabs = sppop2[which(sppop2$year == max(sppop2$year)),]

pmain = ggplot(data = totp,aes(x = year,y = N_med))+
  geom_ribbon(aes(x = year,ymin = N_lci,ymax = N_uci),fill = "#31688EFF",alpha = 0.2)+
  geom_line(color = "#440154FF", linewidth = 1)+
  labs(x = "Year",
       y = "No. of avian nectarivores (Millions)")+
  scale_x_continuous(limits = c(1970, 2020), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal()+
  theme(legend.position = "none", panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 25, margin = margin(t = 5)),
        axis.text.y = element_text(color = "black", size = 25, margin = margin(r = 5)),
        axis.title = element_text(size = 25),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
  )

pdf(paste0("output/Continental level avian pollinator population trajectory.pdf"),
    width = 10,
    height = 7)
print(pmain)
dev.off()
png("output/Continental level avian pollinator population trajectory.png", width = 3000, height = 2500, res = 300)
print(pmain)
dev.off()


# rescaling raw index and smoothed index to ahundance and adding raw ahundance points to fitted line
lost.st = lost.s[rev(order(lost.s$Loss_med)),] #sorting species from largest decrease to largest increase
spord = unique(lost.st$species)
sppop2$spsort = factor(sppop2$species,levels = spord,ordered = T)

sppop2 = sppop2[order(sppop2$spsort,sppop2$year),]

for(cl in c("lci","uci","med","lqrt","uqrt","rescindex","rescindex.raw")){
  sppop2[,cl] = sppop2[,cl]/1e6
}

spinlabs = sppop2[which(sppop2$year == 1970),]
spinlabs = merge(spinlabs,lost.s[,c("species","Loss_med","Loss_lci","Loss_uci","decline")],by = "species")
spinlabs = spinlabs[order(spinlabs$spsort),]

decs = which(spinlabs$decline)
gns = which(spinlabs$decline == F)
spinlabs[decs,"labs"] = paste0(signif(-1*(spinlabs[decs,"Loss_med"])/1e6,2),"M "," [",signif(-1*(spinlabs[decs,"Loss_uci"])/1e6,2)," : ",signif(-1*(spinlabs[decs,"Loss_lci"])/1e6,2),"]")

spinlabs[gns,"labs"] = paste0("+",signif(-1*(spinlabs[gns,"Loss_med"])/1e6,2),"M "," [",signif(-1*(spinlabs[gns,"Loss_uci"])/1e6,2)," : ",signif(-1*(spinlabs[gns,"Loss_lci"])/1e6,2),"]")

rwsnodat = which(is.na(sppop2$rescindex))

pdf(paste0("output/Predicted change in avian pollinator populations by species.pdf"))
for(jj in 1:ceiling(nspecies/9)){
  pmain = ggplot(data = sppop2,aes(x = year,y = med))+
    geom_ribbon(data = sppop2,aes(x = year,ymin = lci,ymax = uci, fill = Breeding.Biome),alpha = 0.2)+
    geom_line(data = sppop2,aes(x = year,y = med))+
    scale_fill_manual(values = fig_palette) +
    labs(x = "",y = "No. of avian nectarivores (Millions)")+
    theme_minimal()+
    geom_text(data = spinlabs,aes(x = 1990,y = Inf,label = labs),vjust = 2,size = 5, inherit.aes = FALSE)+
    theme(legend.position = "none",
          
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          axis.line.x = element_line(color = "black", linewidth = 0.6),
          axis.line.y = element_line(color = "black", linewidth = 0.6),
          
          axis.title.y = element_text(size = 18, margin = margin(r = 10)), 
          
          axis.text = element_text(color = "black", size = 15),
          axis.text.x = element_text(margin = margin(t = 15)),
          
          strip.text = element_text(size = 15))+
    
    facet_wrap_paginate(~spsort,ncol = 3,nrow = 3,scales = "free_y",page = jj)
  print(pmain)
  
  ggsave(
    filename = paste0("output/Predicted change in avian pollinator populations by species page ", jj, ".png"), 
    plot = pmain, 
    width = 10, 
    height = 10, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}
dev.off()

write.csv(indices,"output/Avian pollinator original data.csv")

#Above here is the part of species*year abundance estimation

#-------------------------------------------------------------------
#--------------------------------------------------------------------
#-----------------------------------------------------------------------

#extract posterior distribution of N of each species in each year

N_mat <- MCMCchains(jagsMod$samples, params = "N")

N_dt <- as.data.table(N_mat)
N_dt[, draw := .I]
N_long <- melt(N_dt, id.vars = "draw", variable.name = "node", value.name = "N_draw")

idx <- str_match(N_long$node, "^N\\[(\\d+),(\\d+)\\]$")
N_long[, species_idx := as.integer(idx[,2])]
N_long[, year_idx    := as.integer(idx[,3])]
N_long[, node := NULL]

species_vec <- splist$spfactor #Species in alphabetical order
years_vec <- 1970:2017

N_long[, species := species_vec[species_idx]]
N_long[, year    := years_vec[year_idx]]

# # Convert mean/sd on original (positive) scale to lognormal parameters (mu, sigma) on log scale.
# lognorm_mu_sigma <- function(mean, sd) {
#   if (sd <= 0 || is.na(sd)) stop("Lognormal sd must be > 0")
#   s2 <- log(1 + (sd^2)/(mean^2))
#   list(mu = log(mean) - 0.5*s2, sigma = sqrt(s2))
# }

# Summarize a numeric vector with common statistics (returns a named numeric vector)
summarize_draws <- function(x) {
  list(mean = mean(x), sd = sd(x),
       med = as.numeric(quantile(x, 0.5)),
       lci = as.numeric(quantile(x, 0.025)),
       uci = as.numeric(quantile(x, 0.975)),
       lqrt = as.numeric(quantile(x, 0.25)),
       uqrt = as.numeric(quantile(x, 0.75)))
}


#-------------------------species FVR estimate--------------------------------------------
pol_param <- pollinators

pol_param <- merge(
  x = pol_param,                                         
  y = splist[, c("species", "Breeding.Biome", "Winter.Biome")],
  by.x = "species",                                  
  by.y = "species",                                      
  all.x = TRUE                                          
)

pol_param <- pol_param[order(pol_param$species),]

setDT(pol_param)

FVR_list <- vector("list", nrow(pol_param))

N_FVR_DRAWS <- nrow(N_mat) # keep the number of draws for population size and FVR consistent  

mass.sd <- 0.152 

set.seed(3451)

for (i in seq_len(nrow(pol_param))) {
  r <- pol_param[i,]
  
  #assuming body mass follows truncated normal dist with SD being 15% of mean and lower bound being 1g
  ms <- TruncatedNormal::rtnorm(N_FVR_DRAWS, r$mass_mean, mass.sd*r$mass_mean, lb=1, ub=Inf) 
  
  # assuming predicted FMR follows log-normal dist
  if (r$FMR_pred_a <= 0) stop("FMR_pred_a must be > 0")
  log10FMR_mean <- log10(r$FMR_pred_a) + r$FMR_pred_b * log10(ms)  #predicted FMR treated as median
  lnFMR_mean <- log10FMR_mean*log(10)
  log10FMR_sd <- r$c * sqrt(pmax(0, r$d + r$e * (log10(ms) - r$mean.log.x)^2)) / 1.96  #log-scale SD derived from Nagy’s 95% prediction interval
  lnFMR_sd <- log10FMR_sd*log(10)
  
  lnFMR <- rnorm(N_FVR_DRAWS, lnFMR_mean, lnFMR_sd)
  FMR  <- exp(lnFMR) #kJ/day
  
  # assuming FN and FNE follow uniform dist
  FN  <- runif(N_FVR_DRAWS, r$FN_min,  r$FN_max) #fraction of nectar in diet
  FNE <- runif(N_FVR_DRAWS, r$FNE_min, r$FNE_max) #fraction of nectar extracted from flowers
  
  # assuming NV and NC follow log-normal dist
  NV <- rlnorm(N_FVR_DRAWS, r$NV_mean_ln, r$NV_sd_ln)  # nectar volume in uL
  NC <- rlnorm(N_FVR_DRAWS, r$NC_mean_ln, r$NC_sd_ln)   # nectar concentration in g/ml
  
  ul2ml <- 0.001 # converting ul to ml
  EC_sucorse <- 16.74 # energy content of sucrose in kJ/g
  
  EperF <- NV * ul2ml * EC_sucorse * FNE * NC # energy provided by each flower in kJ/flower
  FVR   <- FMR * FN / EperF     # visits per day in flowers/day
  
  FVR_list[[i]] <- data.table(
    species        = r$species,                 # aligned species name
    draw_fvr       = seq_len(N_FVR_DRAWS),      # local draw index
    FVR            = FVR,
    Service_days   = r$Service_days,            # days stayed in North America
    Breeding.Biome = r$Breeding.Biome,          # carried for grouping (if present)
    Order          = r$Order,
    specialization_lvl    = r$specialization_lvl,
    Habitat          = r$Habitat
  )
}
FVR_dt <- rbindlist(FVR_list)

#----------------------------------------------------------------
# Calculating time series of annual flower visits of each species
# FVR*Population abundance*Service days

K_pop <- max(N_long$draw)

FVR_dt[, draw := ((draw_fvr - 1L) %% K_pop) + 1L]

setkey(FVR_dt, species, draw)
setkey(N_long, species, draw)

spyr <- merge(
  N_long[, .(species, year, draw, N_draw)],
  FVR_dt[, .(species, draw, FVR, Service_days)],
  by = c("species","draw"),
  allow.cartesian = TRUE
)


# -------------------------- Species-year summaries -------------------

spyr[, annual.FV := FVR * Service_days * N_draw] # visits per year

# Summaries for both N_draw and annual.FV
spyr_summ <- spyr[, {
  cN   <- summarize_draws(N_draw) #summarize the distribution of population size
  cAFV <- summarize_draws(annual.FV) #summarize the distribution of AFV
  out  <- c(setNames(cN,   paste0("N_",   names(cN))),
            setNames(cAFV, paste0("AFV_", names(cAFV))))
  out
}, by = .(species, year)]

# write.csv(spyr_summ,"output/Species level Population and AFV trajectories.csv")

# -------------------------- Aggregations by levels and Summarization ---------------------------------
# Continental level: sum across all species (year x draw)
cont <- spyr[, .(N_draw = sum(N_draw), annual.FV = sum(annual.FV)), by = .(year, draw)]

cont_summ <- cont[, {
  cN   <- summarize_draws(N_draw)
  cAFV <- summarize_draws(annual.FV)
  out  <- c(setNames(cN,   paste0("N_",   names(cN))),
            setNames(cAFV, paste0("AFV_", names(cAFV))))
  out
}, by = year]

fwrite(cont_summ,"output/Continental level N and AFV trajectories.csv")


#continental AFV trajectory
for(j in c("AFV_med","AFV_lci","AFV_uci")){
  cont_summ[, (j) := get(j) / 1e9]
}

pmain = ggplot(data = cont_summ,aes(x = year,y = AFV_med))+
  geom_ribbon(aes(x = year,ymin = AFV_lci,ymax = AFV_uci),fill = "#31688EFF",alpha = 0.2)+
  geom_line(color = "#440154FF", linewidth = 1)+
  labs(x = "Year",
       y = "AFV (Billions)")+
  scale_x_continuous(limits = c(1970, 2020), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal()+
  theme(legend.position = "none", panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 25, margin = margin(t = 5)),
        axis.text.y = element_text(color = "black", size = 25, margin = margin(r = 5)),
        axis.title = element_text(size = 25),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
  )

pdf(paste0("output/Continental level AFV trajectory.pdf"),
    width = 10,
    height = 7)
print(pmain)
dev.off()
png("output/Continental level AFV trajectory.png", width = 3000, height = 2500, res = 300)
print(pmain)
dev.off()


# Build a species-to-level key for each level
levels <- c("species","Breeding.Biome","Order","specialization_lvl", "Habitat")
keys <- lapply(levels, function(col) pol_param[, .(species, level = get(col))])
names(keys) <- levels

# aggregate group-level draws (level x year x draw)
group_draws <- list()
for (nm in names(keys)) {
  key <- keys[[nm]]
  if (is.null(key)) next
  spyr_g <- merge(spyr, key, by = "species", all.x = TRUE)
  grp    <- spyr_g[, .(N_draw = sum(N_draw), annual.FV = sum(annual.FV)),
                   by = .(level, year, draw)]
  group_draws[[nm]] <- grp
}

# Summarize group-level draws (level x year) and export as csv
summarized_groups <- function(dt, group_name) {
  dt_summary <- dt[, {
    cN   <- summarize_draws(N_draw)
    cAFV <- summarize_draws(annual.FV)
    out  <- c(setNames(cN,   paste0("N_",   names(cN))),
              setNames(cAFV, paste0("AFV_", names(cAFV))))
    out
  }, by = .(level, year)] 
  
  file_name <- paste0("output/", group_name, " level N and AFV trajectories.csv")
  fwrite(dt_summary,file_name,row.names = FALSE)
  return(dt_summary)
}

group_summ <- mapply(summarized_groups, 
                     dt = group_draws, 
                     group_name = names(group_draws),
                     SIMPLIFY = FALSE)


# -------------------------- Relative-to-first-year trajectories --------------------
Change_trajectory <- function(dt, value_col, by_cols) {
  base <- dt[year == base.yr, .(base_val = get(value_col)), by = c(by_cols, "draw")]
  # change = current - base
  x <- merge(dt, base, by = c(by_cols, "draw"))
  x[, rel_change := get(value_col) - base_val]
  # Summarize per group x year
  summ <- x[, as.list(summarize_draws(rel_change)), by = c(by_cols, "year")]
  setnames(summ, old = names(summ)[!(names(summ) %in% c(by_cols, "year"))],
           new = c("mean","sd","med","lci","uci","lqrt","uqrt"))
  summ
}



# Continental change trajectories (population and annual.FV), relative to the first year
cont_chg_N   <- Change_trajectory(cont,      value_col = "N_draw",    by_cols = character(0))
cont_chg_AFV <- Change_trajectory(cont,      value_col = "annual.FV", by_cols = character(0))

cont_N <- copy(cont_chg_N)
cont_AFV <- copy(cont_chg_AFV)

cont_N_export <- copy(cont_chg_N)
cont_AFV_export <- copy(cont_chg_AFV)

name_cols <- setdiff(names(cont_N_export), c("year"))

setnames(cont_N_export,   name_cols, paste0("N_", name_cols)) 
setnames(cont_AFV_export, name_cols, paste0("AFV_", name_cols)) 

combined_cont <- merge(cont_N_export, cont_AFV_export, by = c("year"))

fwrite(combined_cont,"output/Continental level N and AFV loss trajectories.csv")

for(qt in c("lci","uci","med","lqrt","uqrt")){
  cont_N[, (qt) := get(qt) / 1e6]
}

pmain = ggplot(data = cont_N, aes(x = year, y = med)) +
  geom_ribbon(aes(x = year, ymin = lci, ymax = uci), alpha = 0.2, fill = "#31688EFF") +
  geom_line(color = "#440154FF", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  labs(
    x = "Year",
    y = "Change in no. of avian nectarivores (Millions)") + 
  scale_x_continuous(limits = c(1970, 2020), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = seq(-60, 0, by = 10), expand = expansion(mult = c(0.05, 0.1))) +
  theme_minimal()+
  theme(legend.position = "none", panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 25, margin = margin(t = 5)),
        axis.text.y = element_text(color = "black", size = 25, margin = margin(r = 5)),
        axis.title = element_text(size = 25),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
  )

pdf(paste0("output/Continental level change in avian pollinator population.pdf"),
    width = 10,
    height = 7)
print(pmain)
dev.off()
png("output/Continental level change in avian pollinator population.png", width = 3000, height = 2500, res = 300)
print(pmain)
dev.off()

for(qt in c("lci","uci","med","lqrt","uqrt")){
  cont_AFV[, (qt) := get(qt) / 1e9]
}

pmain = ggplot(data = cont_AFV, aes(x = year, y = med)) +
  geom_ribbon(aes(x = year, ymin = lci, ymax = uci), alpha = 0.2, fill = "#31688EFF") +
  geom_line(color = "#440154FF", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  labs(
    x = "Year",
    y = "Change in AFV (Billions)") + 
  scale_x_continuous(limits = c(1970, 2020), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = pretty, expand = c(0, 0.5)) +
  theme_minimal()+
  theme(legend.position = "none", panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 25, margin = margin(t = 5)),
        axis.text.y = element_text(color = "black", size = 25, margin = margin(r = 5)),
        axis.title = element_text(size = 25),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
  )

pdf(paste0("output/Continental level change in AFV.pdf"),
    width = 10,
    height = 7)
print(pmain)
dev.off()
png("output/Continental level change in AFV.png", width = 3000, height = 2500, res = 300)
print(pmain)
dev.off()

#----------------------------------

# Group-level change trajectories
skip_groups <- c("species", "Order")
group_chg <- list()

for (nm in names(group_draws)) {
  if (nm %in% skip_groups) {
    dt <- group_draws[[nm]]
    chg_N   <- Change_trajectory(dt, value_col = "N_draw",    by_cols = "level")
    chg_AFV <- Change_trajectory(dt, value_col = "annual.FV", by_cols = "level")
    group_chg[[nm]] <- list(N = chg_N, AFV = chg_AFV)
    next
  }
  cat("Processing group:", nm, "\n")
  
  dt <- group_draws[[nm]]
  chg_N   <- Change_trajectory(dt, value_col = "N_draw",    by_cols = "level")
  chg_AFV <- Change_trajectory(dt, value_col = "annual.FV", by_cols = "level")
  group_chg[[nm]] <- list(N = chg_N, AFV = chg_AFV)
  
  export_N   <- copy(chg_N)
  export_AFV <- copy(chg_AFV)
  
  stat_cols <- setdiff(names(export_N), c("level", "year"))
  
  setnames(export_N,   stat_cols, paste0("N_", stat_cols)) # e.g., N_med
  setnames(export_AFV, stat_cols, paste0("AFV_", stat_cols)) # e.g., AFV_med
  
  combined_dt <- merge(export_N, export_AFV, by = c("level", "year"))
  
  file_name_combined <- paste0("output/", nm, " level N and AFV loss trajectories.csv")
  fwrite(combined_dt,file_name_combined,row.names = FALSE)
  
  #-------------------------------------------------
  chg_AFV_level <- copy(chg_AFV)
  for(qt1 in c("lci","uci","med","lqrt","uqrt")){
    chg_AFV_level[, (qt1) := get(qt1) / 1e9]
  }
  
  label_dt1 <- chg_AFV_level[year == 2010]
  
  pmain = ggplot(data = chg_AFV_level,aes(x = year,y = med))+
    geom_ribbon(aes(x = year,ymin = lci,ymax = uci,group = level,fill = level),alpha = 0.2)+
    geom_line(aes(colour = level))+
    labs(x = "Year",y = "Change in AFV (Billions)")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by=10), expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(breaks = pretty, expand = c(0, 0.5)) +
    scale_fill_manual(values = fig_palette) + 
    scale_color_manual(values = fig_palette)+
    theme_minimal()+
    theme(legend.position = "none", panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor = element_blank(), 
          axis.line.x = element_line(color = "black", linewidth = 0.6),
          axis.line.y = element_line(color = "black", linewidth = 0.6),
          axis.text.x = element_text(color = "black", size = 25, margin = margin(t = 5)),
          axis.text.y = element_text(color = "black", size = 25, margin = margin(r = 5)),
          axis.title = element_text(size = 25),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
    )
  
  pdf(paste0("output/", nm, " level AFV change trajectories.pdf"), width = 10, height = 7)
  print(pmain)
  dev.off()
  
  png(paste0("output/", nm, " level AFV change trajectories.png"), width = 3000, height = 2500, res = 300)
  print(pmain)
  dev.off()
  
}


# -------------------------- Plotting AFV trajectories by species ------------------------

spafvchg <- group_chg[["species"]][["AFV"]]
colnames(spafvchg)[colnames(spafvchg) == "level"] <- "species"
spafvchg_2017 <- spafvchg[year==2017]
spafvchg_2017$AFV_decline = F
spafvchg_2017[which(spafvchg_2017$med < 0),"AFV_decline"] = T

spafv <- copy(spyr_summ)
spafv <- merge(spafv, splist3[,c("species", "Breeding.Biome")], by=c("species"), all.x = TRUE)

spafvchg_2017 <- spafvchg_2017[order(spafvchg_2017$med),] #sorting species from largest decrease to largest increase
spafvchg_2017_ord = unique(spafvchg_2017$species)

spafvchg_2017$AFV_sort = factor(spafvchg_2017$species,levels = spafvchg_2017_ord,ordered = T)
spafv$AFV_sort = factor(spafv$species,levels = spafvchg_2017_ord,ordered = T)

spafv = spafv[order(spafv$AFV_sort,spafv$year),]

for(j in c("AFV_med","AFV_lci","AFV_uci")){
  spafv[, (j) := get(j) / 1e9]
}

decs_afv = which(spafvchg_2017$AFV_decline)
gns_afv = which(spafvchg_2017$AFV_decline == F)

spafvchg_2017[decs_afv, AFV_labs := paste0(
  signif(med/1e9, 2), "B ", 
  " [", signif(lci/1e9, 2), " : ", signif(uci/1e9, 2), "]"
)]

spafvchg_2017[gns_afv, AFV_labs := paste0(
  "+", signif(med/1e9, 2), "B ", 
  " [", signif(lci/1e9, 2), " : ", signif(uci/1e9, 2), "]"
)]

pdf(paste0("output/Predicted change in AFV by species.pdf"))
for(jj in 1:ceiling(nspecies/9)){
  pmain = ggplot(data = spafv,aes(x = year,y = AFV_med))+
    geom_ribbon(data = spafv,aes(x = year,ymin = AFV_lci,ymax = AFV_uci, fill = Breeding.Biome),alpha = 0.2)+
    geom_line(data = spafv,aes(x = year,y = AFV_med))+
    scale_fill_manual(values = fig_palette) +
    labs(x = "",y = "AFV (Billions)")+
    theme_minimal()+
    geom_text(data = spafvchg_2017,aes(x = 1990,y = Inf,label = AFV_labs),vjust = 2,size = 5, inherit.aes = FALSE)+
    theme(legend.position = "none",
          
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          axis.line.x = element_line(color = "black", linewidth = 0.6),
          axis.line.y = element_line(color = "black", linewidth = 0.6),
          
          axis.title.y = element_text(size = 18, margin = margin(r = 10)), 
          
          axis.text = element_text(color = "black", size = 15),
          axis.text.x = element_text(margin = margin(t = 15)),
          
          strip.text = element_text(size = 15))+
    
    facet_wrap_paginate(~AFV_sort,ncol = 3,nrow = 3,scales = "free_y",page = jj)
  print(pmain)
  
  ggsave(
    filename = paste0("output/Predicted change in AFV by species page ", jj, ".png"), 
    plot = pmain, 
    width = 10, 
    height = 10, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}
dev.off()


# -------------------------- AFV loss summary table --------------------------

AFV_prop_loss <- function(dt, value_col, by_cols) {
  # Extract base year values (1970)
  base <- dt[year == base.yr, .(base_val = get(value_col)), by = c(by_cols, "draw")]
  # Join current data with base data
  x <- merge(dt[year == end.yr], base, by = c(by_cols, "draw"))
  # Calculate proportional loss: (2017 - 1970) / 1970
  x[, p_loss := (get(value_col) - base_val) / base_val]
  # Summarize the distribution of proportional loss
  summ <- x[, as.list(summarize_draws(p_loss)), by = c(by_cols)]
  # Rename columns
  setnames(summ, 
           old = c("mean", "sd", "med", "lci", "uci", "lqrt", "uqrt"),
           new = c("plostmean", "plostsd", "plostmed", "plostlci", "plostuci", "plostlqrt", "plostuqrt"))
  return(summ)
}

# Calculate the difference: AFV_2017 - AFV_1970 for every species
sp_chg_AFV <- Change_trajectory(spyr, value_col = "annual.FV", by_cols = "species")

# Extract the final year results to identify species-level "losers"
lost.afv.s <- sp_chg_AFV[year == 2017, .(species, med)]
setnames(lost.afv.s, "med", "AFV_loss_med")

# Define decline as a negative net change
lost.afv.s[, decline := AFV_loss_med < 0]

# Merge with metadata for grouping
lost.afv.s <- merge(lost.afv.s, 
                    pol_param[, .(species, Breeding.Biome, Order, specialization_lvl, Habitat)], 
                    by = "species")

# Define levels to include in the summary (excluding 'species' as it is the base)
sum_levels <- c("Breeding.Biome", "specialization_lvl", "Habitat")

# Initialize with Continental level change (Total AFV change)
cont_AFV_lost <- cont_AFV[year == end.yr]
cont_AFV_plost <- AFV_prop_loss(cont, value_col = "annual.FV", by_cols = character(0))
afv_alll <- data.frame(
  Group    = "Continental",
  nspecies = nspecies,
  lci      = round(cont_AFV_lost$lci, 1),
  med      = round(cont_AFV_lost$med, 1),
  uci      = round(cont_AFV_lost$uci, 1),
  plostlci = round(cont_AFV_plost$plostlci, 2),
  plostmed = round(cont_AFV_plost$plostmed, 2),
  plostuci = round(cont_AFV_plost$plostuci, 2),
  node     = "AFVlost.all" 
)

afv_list <- list(afv_alll)

# Loop through each group level to extract summary stats
for (nm in sum_levels) {
  
  # Get absolute net change (med, lci, uci)
  dt_abs <- copy(group_chg[[nm]]$AFV[year == end.yr])
  
  # Get proportional loss change (plostmed, plostlci, plostuci)
  dt_prop <- AFV_prop_loss(group_draws[[nm]], value_col = "annual.FV", by_cols = "level")
  
  # Merge stats
  dt_lvl <- merge(dt_abs, dt_prop, by = "level")
  
  # Count species within this group
  sp_count <- lost.afv.s[, .(nspp = .N), by = c(nm)]
  setnames(sp_count, nm, "level")
  
  dt_lvl <- merge(dt_lvl, sp_count, by = "level")
  
  df_lvl <- data.frame(
    Group    = dt_lvl$level,
    nspecies = dt_lvl$nspp,
    lci      = round(dt_lvl$lci/1e9,1),
    med      = round(dt_lvl$med/1e9,1),
    uci      = round(dt_lvl$uci/1e9,1),
    plostlci = round(dt_lvl$plostlci,2),
    plostmed = round(dt_lvl$plostmed,2),
    plostuci = round(dt_lvl$plostuci,2),
    # Generate node names
    node     = paste0("AFVlost.", tolower(gsub("\\.", "", nm)), "[", 1:nrow(dt_lvl), "]")
  )
  afv_list[[nm]] <- df_lvl
}

# Combine all rows
afvsout <- do.call(rbind, afv_list)

# calculating the number and proportion of species with declining AFV in each category
afvsout$proportion.species.AFVdecline <- NA
afvsout$n.species.AFVdecline <- NA

for(i in 1:nrow(afvsout)){
  current_group <- afvsout[i, "Group"]
  current_node  <- afvsout[i, "node"]
  
  if(current_group == "Continental"){
    S_decline = lost.afv.s$decline
  } else {
    # Extract level name from node string to find the correct metadata column
    lvl_type <- gsub("AFVlost\\.", "", gsub("\\[.*", "", current_node))
    
    # Mapping table for metadata columns
    col_map <- c(breedingbiome = "Breeding.Biome", order = "Order", 
                 specialization_lvl = "specialization_lvl", habitat = "Habitat")
    meta_col <- col_map[lvl_type]
    
    # Get status for species belonging to this specific group
    S_decline = lost.afv.s[get(meta_col) == current_group, decline]
  }
  
  if(length(S_decline) > 0){
    afvsout[i, "proportion.species.AFVdecline"] <- round(sum(S_decline) / length(S_decline), 1)
    afvsout[i, "n.species.AFVdecline"]          <- sum(S_decline)
  }
}

fwrite(afvsout, "output/Summary-Net change in AFV across the North American avian pollinators.csv", row.names = F)


# ------------------ AFV loss Forest Plot for biomes and habitats -------------------

AFVloss_forest <- afvsout[afvsout$Group != "Continental" & 
                            (grepl("breedingbiome", afvsout$node) | grepl("habitat", afvsout$node)), ]

# Create a 'Level' column to distinguish between Biome and Habitat for facetting
AFVloss_forest$Level <- ifelse(grepl("breedingbiome", AFVloss_forest$node), "Breeding Biome", "Habitat")

forest_plot <- function(data) {
  
  # # Sort data by median loss (plostmed) from highest to lowest
  # data$original_order <- 1:nrow(data)
  #data$Group <- reorder(data$Group, data$plostmed)
  
  ggplot(data, aes(x = plostmed*100, y = Group, color = Group)) +
    # Add a vertical line at 0 (no change)
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.6) +
    geom_pointrange(aes(xmin = plostlci*100, xmax = plostuci*100), size = 1, linewidth = 1) +
    labs(x = "Change in AFV since 1970 (%)",
         y = "") +
    scale_x_continuous(breaks = pretty, expand = expansion(mult = c(0.05, 0.1))) +
    scale_color_manual(values = fig_palette)+
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      panel.grid.major = element_line(color = "gray50", linewidth = 0.2, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(), 
      axis.text.y = element_text(color = "black", size = 25),
      axis.text.x = element_text(color = "black", size = 25),
      axis.title.x = element_text(size = 25, margin = margin(t = 10))
    )
}

biome_data <- subset(AFVloss_forest, Level == "Breeding Biome")
p_biome <- forest_plot(biome_data)

habitat_data <- subset(AFVloss_forest, Level == "Habitat")
p_habitat <- forest_plot(habitat_data)

png("output/Breeding Biome level AFV proportional change forest plot.png", 
    width = 2250, height = 2500, res = 300)
print(p_biome)
dev.off()
pdf("output/Breeding Biome level AFV proportional change forest plot.pdf", 
    width = 5, height = 7)
print(p_biome)
dev.off()

png("output/Habitat level AFV proportional change forest plot.png", 
    width = 2000, height = 2500, res = 300)
print(p_habitat)
dev.off()
pdf("output/Habitat level AFV proportional change forest plot.pdf", 
    width = 4, height = 7)
print(p_habitat)
dev.off()


# ----- Percentages of each species' contribution to total AFV in 1970 and 2017--------

cont_totals <- cont[year %in% c(base.yr, end.yr), .(year, draw, total_AFV = annual.FV)]
sp_sub <- spyr[year %in% c(base.yr, end.yr), .(species, year, draw, sp_AFV = annual.FV)]
contribution_draws <- merge(sp_sub, cont_totals, by = c("year", "draw"))
contribution_draws[, pct_contrib := sp_AFV / total_AFV*100]

# contribution of specialized and generalized nectarivores to ttl AFV in 1970 and 2017
contribution_draws_spcl <- merge(contribution_draws, splist3[, c("species", "specialization_lvl")], by = "species", all.x = TRUE)

contrib_spcl <- contribution_draws_spcl[, .(
  pct_contrib_spcl = sum(pct_contrib, na.rm = TRUE)
), by = .(draw, year, specialization_lvl)]

contrib_summ_spcl <- contrib_spcl[, {
  pct <- summarize_draws(pct_contrib_spcl)
  pct
}, by = .(year, specialization_lvl)]

# each species' contribution to total AFV in 1970 and 2017
contrib_summ <- contribution_draws[, {
  pct <- summarize_draws(pct_contrib)
  pct
}, by = .(species, year)]

top_n <- 20 # only display the top 20 species

contrib_1970 <- contrib_summ[year == 1970]
contrib_1970 <- contrib_1970[order(med)][1:top_n]
contrib_2017 <- contrib_summ[year == 2017]
contrib_2017 <- contrib_2017[order(med)][1:top_n]

png("output/Percentages of species-level contribution to total AFV in 1970 and 2017.png", 
    width = 14, height = 8, units = "in", res = 300, pointsize = 12)

op <- par(mfrow = c(1, 2), mar = c(3, 16, 1, 1), oma = c(2, 0, 0, 0))

bp_70 <- barplot(height = contrib_1970$med, 
                 names.arg = contrib_1970$species,
                 cex.names=1.33, horiz=TRUE, las=1, xlim=c(0,40), cex.axis=1.33,
                 main = "",
                 space=0.15,
                 ylim = c(0,20),
                 col = "lightgrey", 
                 border = "black",
                 ylab = "")

x_ticks_70 <- seq(10, 30, by = 10)

segments(x0 = x_ticks_70, 
         y0 = min(bp_70) - 1.3, 
         x1 = x_ticks_70, 
         y1 = max(bp_70)-3.5, 
         col = "lightgrey", lty = "dotted")

arrows(y0 = bp_70, x0 = contrib_1970$lci, 
       y1 = bp_70, x1 = contrib_1970$uci, 
       angle = 90, code = 3, length = 0.04, col = "black", lwd = 1.5)

text(x = 0, y = 20.5, labels = "(a) 1970", cex = 1.33, xpd=TRUE, adj = 0)

bp_17 <- barplot(height = contrib_2017$med, 
                 names.arg = contrib_2017$species,
                 cex.names=1.33, horiz=TRUE, las=1, xlim=c(0,50), cex.axis=1.33,
                 main = "",
                 space=0.15,
                 ylim = c(0,20),
                 col = "lightgrey", 
                 border = "black",
                 ylab = "")

x_ticks_17 <- seq(10, 40, by = 10)

segments(x0 = x_ticks_17, 
         y0 = min(bp_17) - 1.3, 
         x1 = x_ticks_17, 
         y1 = max(bp_17)-3.5, 
         col = "lightgrey", lty = "dotted")

arrows(y0 = bp_17, x0 = contrib_2017$lci, 
       y1 = bp_17, x1 = contrib_2017$uci, 
       angle = 90, code = 3, length = 0.04, col = "black", lwd = 1.5)

text(x = 0, y = 20.5, labels = "(b) 2017", cex = 1.33, xpd=TRUE, adj = 0)

mtext("Species-level contribution to total AFV in North America (%)", 
      side = 1, line = 0.4, outer = TRUE, cex = 1.5)

par(op)
dev.off()


# ---------------------- Histogram of species-level relative change in AFV 1970 to 2017 ------------------------

sp_70_17 <- dcast(sp_sub, species + draw ~ year, value.var = "sp_AFV")
setnames(sp_70_17, c("species", "draw", "AFV_1970", "AFV_2017"))
sp_70_17 <- merge(sp_70_17, splist3[, c("species", "Breeding.Biome", "Habitat")], by = "species", all.x = TRUE)

sp_70_17[, rel_change := (AFV_2017 - AFV_1970) / AFV_1970 * 100]

sp_rel_summ <- sp_70_17[, {
  rel <- summarize_draws(rel_change)
  rel
}, by = species]

sp_rel_summ[, change_type := ifelse(med < 0, "Decrease", "Increase")]

sp_rel_summ$species <- reorder(sp_rel_summ$species, sp_rel_summ$med)

color_palette <- c("Decrease" = "#E69F00", "Increase" = "#56B4E9")

p_rel <- ggplot(sp_rel_summ, aes(x = med, y = species, fill = change_type)) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid",  linewidth = 0.6) +
  geom_col(width = 0.8, alpha = 0.85) +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.2, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(breaks=pretty) +
  labs(
    x="Change in AFV since 1970 (%)",
    y=NULL
  )+
  theme_minimal()+
  theme(legend.position = "none", 
        panel.grid.major = element_line(color = "gray50", linewidth = 0.2, linetype = "dashed"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(color = "black", size = 14, margin = margin(r = 5)),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10)),
  )

png("output/Histogram of species-level relative change in AFV 1970 to 2017.png", 
    width = 1800, height = 2000, res = 300)
print(p_rel)
dev.off()
pdf("output/Histogram of species-level relative change in AFV 1970 to 2017.pdf", 
    width = 6, height = 8)
print(p_rel)
dev.off()

# --------------- Histogram of species-level absolute change in AFV 1970 to 2017 by biome ------------------------

sp_70_17[, abs_change := AFV_2017 - AFV_1970]

sp_abs_summ <- sp_70_17[, {
  abs <- summarize_draws(abs_change)
  abs
}, by = c("species", "Breeding.Biome", "Habitat")]

sp_abs_summ[, change_type := ifelse(med < 0, "Decrease", "Increase")]

for(qt in c("lci","uci","med")){
  sp_abs_summ[, (qt) := get(qt) / 1e9]
}

sp_abs_summ_biome <- sp_abs_summ[order(Breeding.Biome, med)]
sp_abs_summ_biome$species <- factor(sp_abs_summ_biome$species, levels = sp_abs_summ_biome$species)

color_palette <- c("Decrease" = "#E69F00", "Increase" = "#56B4E9")

p_abs <- ggplot(sp_abs_summ_biome, aes(x = med, y = species, fill = change_type)) +
  
  geom_col(width = 0.8, alpha = 0.85) +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.2, color = "black", linewidth = 0.4) +
  
  geom_text(aes(
    x = ifelse(med < 0, 30, -30), label = species, hjust = ifelse(med < 0, 0, 1)), size = 4.2, color = "black") +
  
  geom_vline(xintercept = 0, color = "black", linetype = "solid",  linewidth = 0.6) +
  scale_fill_manual(values = color_palette) +
  
  facet_grid(Breeding.Biome ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  scale_x_continuous(limits = c(-450, 600), breaks = seq(-400, 600, by = 200)) +
  
  labs(
    x="Absolute change in AFV since 1970 (Billions)",
    y=NULL
  )+
  theme_minimal()+
  theme(legend.position = "none", 
        panel.grid.major = element_line(color = "gray50", linewidth = 0.2, linetype = "dashed"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 14, margin = margin(t = 5)),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10)),
        strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1, size = 16), 
        strip.placement = "outside", 
        panel.spacing = unit(1, "lines") 
  )

png("output/Histogram of species-level absolute change in AFV 1970 to 2017 by biome.png", 
    width = 2200, height = 2200, res = 300)
print(p_abs)
dev.off()
pdf("output/Histogram of species-level absolute change in AFV 1970 to 2017 by biome.pdf", 
    width = 8, height = 8)
print(p_abs)
dev.off()

# --------------- Histogram of species-level absolute change in AFV 1970 to 2017 by habitat ------------------------

sp_abs_summ_hab <- sp_abs_summ[order(Habitat, med)]
sp_abs_summ_hab$species <- factor(sp_abs_summ_hab$species, levels = sp_abs_summ_hab$species)

p_abs <- ggplot(sp_abs_summ_hab, aes(x = med, y = species, fill = change_type)) +
  
  geom_col(width = 0.8, alpha = 0.85) +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.2, color = "black", linewidth = 0.4) +
  
  geom_text(aes(
    x = ifelse(med < 0, 30, -30), label = species, hjust = ifelse(med < 0, 0, 1)), size = 4.2, color = "black") +
  
  geom_vline(xintercept = 0, color = "black", linetype = "solid",  linewidth = 0.6) +
  scale_fill_manual(values = color_palette) +
  
  facet_grid(Habitat ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  scale_x_continuous(limits = c(-450, 600), breaks = seq(-400, 600, by = 200)) +
  
  labs(
    x="Absolute change in AFV since 1970 (Billions)",
    y=NULL
  )+
  theme_minimal()+
  theme(legend.position = "none", 
        panel.grid.major = element_line(color = "gray50", linewidth = 0.2, linetype = "dashed"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.text.x = element_text(color = "black", size = 14, margin = margin(t = 5)),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10)),
        strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1, size = 16), 
        strip.placement = "outside", 
        panel.spacing = unit(1, "lines") 
  )

png("output/Histogram of species-level absolute change in AFV 1970 to 2017 by habitat.png", 
    width = 2000, height = 2200, res = 300)
print(p_abs)
dev.off()
pdf("output/Histogram of species-level absolute change in AFV 1970 to 2017 by habitat.pdf", 
    width = 8, height = 8)
print(p_abs)
dev.off()

# --------------- boxplot of AFV change distribution by habitat in each biome------------------------

grp_70_17 <- sp_70_17[, .(
  grp_AFV_1970 = sum(AFV_1970, na.rm = TRUE),
  grp_AFV_2017 = sum(AFV_2017, na.rm = TRUE)
), by = .(Breeding.Biome, Habitat, draw)]

grp_70_17[, grp_rel_change := (grp_AFV_2017 - grp_AFV_1970) / grp_AFV_1970 * 100] 

grp_rel_summ <- grp_70_17[, {
  rel <- summarize_draws(grp_rel_change)
  rel
}, by = c("Breeding.Biome", "Habitat")]

p_group_box <- ggplot(grp_rel_summ, aes(x = "")) +
  geom_boxplot(
    aes(ymin = lci, lower = lqrt, middle = med, upper = uqrt, ymax = uci),
    stat = "identity",
    width = 0.5,       
    alpha = 1,      
    color = "black",   
    size = 0.4         
  ) +
  facet_grid(rows = vars(Habitat), cols = vars(Breeding.Biome), scales = "fixed") +
  
  labs(
    x = "",
    y = "Relative change in AFV since 1970 (%)"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.6),
    
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 14),
    
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(size = 14, face = "bold"),
    strip.placement = "outside",
    
    legend.position = "none"
  )

png("output/Boxplot of AFV relative change distribution by habitat in each biome.png", width = 1800, height = 2400, res = 300)
print(p_group_box)
dev.off()

pdf("output/Boxplot of AFV relative change distribution by habitat in each biome.pdf", width = 6, height = 8)
print(p_group_box)
dev.off()

# ----- Percentages of each species' contribution to biome-level ttl AFV in 1970 and 2017--------

sp_sub_biome <- merge(sp_sub, pol_param[, .(species, Breeding.Biome)], by = "species")

biome_AFV_ttl <- sp_sub_biome[, .(biome_total = sum(sp_AFV)), by = .(year, draw, Breeding.Biome)]

sp_contrib_biome <- merge(sp_sub_biome, biome_AFV_ttl, by = c("year", "draw", "Breeding.Biome"))
sp_contrib_biome[, pct_contrib := sp_AFV / biome_total*100]

biome_contrib_summ <- sp_contrib_biome[, {
  biome_afv <- summarize_draws(pct_contrib)
  biome_afv
}, by = .(species, year, Breeding.Biome)]

biomes_list <- sort(unique(biome_contrib_summ$Breeding.Biome))
n_biomes <- length(biomes_list)

png("output/Percentages of species-level contribution to biome-level total AFV in 1970 and 2017.png", 
    width = 20, height = 8 * n_biomes, units = "in", res = 300,, pointsize = 12)

mat <- matrix(1:(n_biomes * 3), ncol = 3, byrow = TRUE)
layout(mat, widths = c(1.5, 5, 5))

top_n_biome <- 10

par(oma = c(5, 1, 2, 1), cex = 1.1)

for (i in seq_along(biomes_list)) {
  
  b <- biomes_list[i]
  afv_biome <- biome_contrib_summ[Breeding.Biome == b]
  
  afv_biome70 <- afv_biome[year == 1970][order(med)][tail(seq_len(.N), top_n_biome)]
  afv_biome17 <- afv_biome[year == 2017][order(med)][tail(seq_len(.N), top_n_biome)]
  
  xmax <- max(c(afv_biome70$uci, afv_biome17$uci), na.rm = TRUE) * 1.05
  
  par(mar = c(0, 0, 0, 0)) 
  plot.new()
  text(x = 0.5, y = 0.5, labels = b, cex = 1.8, font = 2, adj = 0.5)
  
  par(mar = c(3, 15, 2, 1), lwd = 3)
  bp70 <- barplot(height = afv_biome70$med, 
                  names.arg = afv_biome70$species,
                  horiz = TRUE, 
                  xlim = c(0, xmax),
                  las = 1, 
                  cex.names = 1.5,
                  cex.axis = 1.5,
                  col = "lightgrey", border = "black")
  
  arrows(y0 = bp70, x0 = afv_biome70$lci, y1 = bp70, x1 = afv_biome70$uci, 
         angle = 90, code = 3, length = 0.08, col = "black", lwd = 2.5)
  
  
  if (i == 1) {
    mtext("(a) 1970", side =3, line =1, cex = 2, xpd=TRUE, adj =0)
  }
  
  par(mar = c(3, 15, 2, 1))
  bp17 <- barplot(height = afv_biome17$med, 
                  names.arg = afv_biome17$species,
                  horiz = TRUE, 
                  xlim = c(0, xmax),
                  las = 1, 
                  cex.names = 1.5,
                  cex.axis = 1.5,
                  col = "lightgrey", border = "black")
  
  arrows(y0 = bp17, x0 = afv_biome17$lci, y1 = bp17, x1 = afv_biome17$uci, 
         angle = 90, code = 3, length = 0.08, col = "black", lwd = 2.5)
  
  if (i == 1) {
    mtext("(b) 2017", side =3, line =1, cex = 2, xpd=TRUE, adj =0)
  }
}

mtext("Species-level contribution to total AFV within each breeding biome (%)", side = 1, line = 1.2, outer = TRUE, cex = 2)

dev.off()


#------------------------Source of uncertainty analysis for AFV 1970 and 2017--------------------------------------------------

summarize_draws_90 <- function(x) {
  list(mean = mean(x), sd = sd(x),
       med = as.numeric(quantile(x, 0.5)),
       lci = as.numeric(quantile(x, 0.05)),
       uci = as.numeric(quantile(x, 0.95)))
}

Change_trajectory_90 <- function(dt, value_col, by_cols) {
  base <- dt[year == base.yr, .(base_val = get(value_col)), by = c(by_cols, "draw")]
  # change = current - base
  x <- merge(dt, base, by = c(by_cols, "draw"))
  x[, rel_change := get(value_col) - base_val]
  # Summarize per group x year
  summ <- x[, as.list(summarize_draws_90(rel_change)), by = c(by_cols, "year")]
  setnames(summ, old = names(summ)[!(names(summ) %in% c(by_cols, "year"))],
           new = c("mean","sd","med","lci","uci"))
  summ
}

cont_summ_90 <- cont[, {
  cN   <- summarize_draws_90(N_draw)
  cAFV <- summarize_draws_90(annual.FV)
  out  <- c(setNames(cN,   paste0("N_",   names(cN))),
            setNames(cAFV, paste0("AFV_", names(cAFV))))
  out
}, by = year]

cont_chg_AFV_90 <- Change_trajectory_90(cont,      value_col = "annual.FV", by_cols = character(0))

source("models/uncertainty source analysis-body mass only.txt")

source("models/uncertainty source analysis-FMR equation only.txt")

source("models/uncertainty source analysis-FN only.txt")

source("models/uncertainty source analysis-FNE only.txt")

source("models/uncertainty source analysis-NV only.txt")

source("models/uncertainty source analysis-NC only.txt")

source("models/uncertainty source analysis-N only.txt")


### plot the contribution of each variable to total uncertainty in AFV change from 1970 to 2017
uncertainty_chg_dt <- data.table(
  
  Source = c("FMR equation", "% of nectar in birds' diet", "% of nectar extracted by birds", 
             "Body Mass", "Nectar concentration", "Nectar volume", "Population size estimates"),
  
  Raw_Chg_Prop = c(
    fmr_chg_prop, fn_chg_prop, fne_chg_prop, 
    ms_chg_prop, nc_chg_prop, nv_chg_prop, n_chg_prop
  )
)

uncertainty_chg_dt[, Normalized_Chg_Prop := Raw_Chg_Prop / sum(Raw_Chg_Prop)*100]

uncertainty_chg_dt[, y_label := "Total AFV Change"]

plot_chg_dt <- copy(uncertainty_chg_dt)

order_source_chg <- plot_chg_dt[order(-Normalized_Chg_Prop)]$Source
plot_chg_dt$Source <- factor(plot_chg_dt$Source, levels = order_source_chg)

plot_chg_uncertainty <- ggplot(plot_chg_dt, aes(y = y_label, x = Normalized_Chg_Prop, fill = Source)) +
  
  geom_col(width = 0.5) +
  geom_text(aes(label = ifelse(Normalized_Chg_Prop > 3, paste0(round(Normalized_Chg_Prop, 1), "%"), "")), 
            position = position_stack(vjust = 0.5), size = 2.5, color = "white", fontface = "bold") +
  
  scale_fill_viridis_d(option = "viridis") +
  
  scale_y_discrete(expand = c(0.1, 0.3)) +
  
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), 
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     expand = c(0, 0)) +
  
  labs(
    title = NULL, x = NULL, y = NULL, fill = NULL
  )+
  
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.length.x = unit(0.2, "cm"),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    axis.text.x = element_text(color = "black", size = 8, vjust = 1),
    
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(t = 0, r = 20, b = 0, l = 20, unit = "pt")
  )

png("output/Source of uncertainty in total AFV change from 1970 to 2017.png", 
    width = 1800, height = 500, res = 300)
print(plot_chg_uncertainty)
dev.off()
pdf("output/Source of uncertainty in total AFV change from 1970 to 2017.pdf", 
    width = 9, height = 2.5)
print(plot_chg_uncertainty)
dev.off()




