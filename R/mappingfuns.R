# Map user space to inches space:
usr2inX <- function(x)
{
  usr <- c(-1,1,-1,1)
  pin <- par("din")
  (x-usr[1])/(usr[2]-usr[1]) * pin[1]
}

usr2inY <- function(x)
{
  usr <- c(-1,1,-1,1)
  pin <- par("din")
  (x-usr[3])/(usr[4]-usr[3]) * pin[2]
}

# Same but about origin (for atan2):
usr2inX2 <- function(x)
{
  usr <- c(-1,1,-1,1)
  pin <- par("din")
  x/(usr[2]-usr[1]) * pin[1]
}

usr2inY2 <- function(x)
{
  usr <- c(-1,1,-1,1)
  pin <- par("din")
  x/(usr[4]-usr[3]) * pin[2]
}
atan2usr2in <- function(x,y) atan2(usr2inX2(x),usr2inY2(y))%%(2*pi)

# Map inches space to user space:
in2usrX <- function(x)
{
  usr <- c(-1,1,-1,1)
  pin <- par("din")
  usr[1] + x/pin[1] * (usr[2] - usr[1])
}

in2usrY <- function(x)
{
  usr <- c(-1,1,-1,1)
  pin <- par("din")
  usr[3] + x/pin[2] * (usr[4] - usr[3])
}