# function calcPMV from comf-package which has thermal load L as output
calcPMVL <- 
function (ta, tr, vel, rh, clo = 0.5, met = 1, wme = 0, basMet = 58.15) 
{
    m <- met * basMet
    w <- wme * basMet
    mw <- m - w
    icl <- 0.155 * clo
    pa <- rh * 10 * exp(16.6536 - (4030.183/(ta + 235)))
    if (icl <= 0.078) {
        fcl <- 1 + 1.29 * icl
    }
    else {
        fcl <- 1.05 + 0.645 * icl
    }
    fcic <- icl * fcl
    p2 <- fcic * 3.96
    p3 <- fcic * 100
    tra <- tr + 273
    taa <- ta + 273
    p1 <- fcic * taa
    p4 <- 308.7 - 0.028 * mw + p2 * (tra/100)^4
    tclA <- taa + (35.5 - ta)/(3.5 * (6.45 * icl + 0.1))
    xn <- tclA/100
    xf <- xn
    hcf <- 12.1 * (vel)^0.5
    noi <- 0
    eps <- 0.00015
    while (noi < 150) {
        xf <- (xf + xn)/2
        hcn <- 2.38 * abs(100 * xf - taa)^0.25
        if (hcf > hcn) {
            hc <- hcf
        }
        else {
            hc <- hcn
        }
        xn <- (p4 + p1 * hc - p2 * xf^4)/(100 + p3 * hc)
        noi <- noi + 1
        if (noi > 1 & abs(xn - xf) <= eps) {
            break
        }
    }
    tcl <- 100 * xn - 273
    pm1 <- 3.96 * fcl * (xn^4 - (tra/100)^4)
    pm2 <- fcl * hc * (tcl - ta)
    pm3 <- 0.303 * exp(-0.036 * m) + 0.028
    if (mw > basMet) {
        pm4 <- 0.42 * (mw - basMet)
    }
    else {
        pm4 <- 0
    }
    pm5 <- 3.05 * 0.001 * (5733 - 6.99 * mw - pa)
    pm6 <- 1.7 * 1e-05 * m * (5867 - pa) + 0.0014 * m * (34 - 
        ta)
    pmv <- pm3 * (mw - pm5 - pm4 - pm6 - pm1 - pm2)
	Lraw <- (mw - pm5 - pm4 - pm6 - pm1 - pm2)
    ppd <- 100 - 95 * exp(-0.03353 * pmv^4 - 0.2179 * pmv^2)
    data.frame(pmv, ppd, Lraw)
}
