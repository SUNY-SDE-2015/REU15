coralDelayModel <- function( y, parameters) {
    M <- y[1]
    C <- y[2]
    a <- parameters[1]
    g <- parameters[2]
    gamma <- parameters[3]
    r <- parameters[4]
    d <- parameters[5]
    dy[1] <- a*M*C - (g*M)/(1-C) + gamma*M - gamma*M*M - gamma*M*C
    dy[2] <- r*C - r*M*C - r*C*C - d*C - a*M*C
    list(dy)
}