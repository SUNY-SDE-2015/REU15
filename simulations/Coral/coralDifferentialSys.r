coralDelayModel <- function(t, y, parameters) {
    M <- y[1]
    C <- y[2]
    a <- parameters[1]
    g <- parameters[2]
    tau <- parameters[3]
    gamma <- parameters[4]
    r <- parameters[5]
    d <- parameters[6]
    dy[1] <- a*M*C - (g*M*(t-tau))/(1-C(t-tau)) + gamma*M - gamma*M*M - gamma*M*C
    dy[2] <- r*C - r*M*C - r*C*C - d*C - a*M*C
    list(dy)
}