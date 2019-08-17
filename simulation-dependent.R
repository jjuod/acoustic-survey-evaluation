#!/usr/bin/Rscript

# USAGE:
#./simulations-singlerun.R D g0 sigma trapnum mu

library(ascr)
args = commandArgs(TRUE)
D = as.numeric(args[1])
g0 = as.numeric(args[2])
sigma = as.numeric(args[3])
trapnum = as.numeric(args[4])
mu = as.numeric(args[5])
print("Simulating DEPENDENT call locations with parameters:")
print(paste("D:", D, "g0:", g0, "sigma:", sigma, "trapnum:", trapnum, "mu:", mu))

# create various options of traps
if(trapnum==3){
    traps = as.matrix(data.frame(x=rep(c(-100, 0, 100), each=3), y=rep(c(-100, 0, 100), 3)))
    mask = create.mask(traps, buffer=200)
} else if(trapnum==4){
    traps = as.matrix(data.frame(x=rep(c(-100, 0, 100, 200), each=4), y=rep(c(-100, 0, 100, 200), 4)))
    mask = create.mask(traps, buffer=200)
} else if(trapnum==5){
    traps = as.matrix(data.frame(x=rep(c(-200, -100, 0, 100, 200), each=5), y=rep(c(-200, -100, 0, 100, 200), 5)))
    mask = create.mask(traps, buffer=200)
}

# for a number of iterations, draw samples from this population,
# get the resulting sampling distribution of the parameters
niter = 100
out = data.frame()
for(i in 1:niter){
    print(paste("working on iter", i))
    capt = sim.capt(traps=traps, mask=mask, pars=list(D=D, g0=g0, sigma=sigma), cue.rates = mu, survey.length = 1)
    test.fit = fit.ascr(capt=capt, traps=traps, mask=mask, detfn="hn", survey.length=1, hess=F, cue.rates = mu)
    res = c(i=i, trueD=D, trueg0=g0, trapnum=trapnum, truesigma=sigma, coef(test.fit))
    out = rbind(out, res)
    names(out) = names(res)
}

write.table(out, paste0("resdep-", D, "-", g0, "-", sigma, "-", trapnum, "-", mu, ".csv"), append=F, quote=F, row.names=F, col.names=T)
