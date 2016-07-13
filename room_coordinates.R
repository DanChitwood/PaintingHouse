#Read in ggplot2 to make graphs and mgcv to find intersection of polygon and points using "in.out"

library(ggplot2)
library(mgcv)

#Read in data for polygon of house

data <- read.table("./intersection.txt", header=TRUE)

#Create a formula for phyllotaxy based on the Fibonacci angle
#The below is modified from http://stackoverflow.com/questions/26141237/how-to-create-phyllotaxis-spirals-with-r, which cites http://algorithmicbotany.org/papers/abop/abop-ch4.pdf

golden.ratio = (sqrt(5) + 1)/2
fibonacci.angle=360/(golden.ratio^2)
c=6.558
num_points=100000
num=rep(0,num_points)
x=rep(0,num_points)
y=rep(0,num_points)
t=rep(0,num_points)

for (n in 1:num_points) {
    r=c*sqrt(n)
    theta=fibonacci.angle*(n)
    num[n]=n
    x[n]=r*cos(theta)
    y[n]=r*sin(theta)
    t[n]=theta
}

plot(x,y,axes=FALSE,ann=FALSE,pch=19,cex=0.01)

#Create a dataframe with the order, x & y coordinates of points, and their angle

spiral <- as.data.frame(cbind(num,x,y,t))

#Determine which points fall within the polygons of the house to be painted
#First, create two matrices, one for the polygons and one for the points

polygons <- as.matrix(data[,6:7])
points <- as.matrix(spiral[,2:3])

#Use the in.out function to find those points within the polygons

inside <- in.out(polygons, points)

#Combine points with other phyllotaxy information

inside_points <- cbind(inside, spiral)

#Keep only points of interest in the polygons

just_inside <- subset(inside_points, inside==TRUE)

#Read in landmark data from 3,319 Passiflora leaves. For more info, see: https://github.com/DanChitwood/PassifloraLeaves

coord <- read.table("./3319_Procrustes_Leaves.txt", header=TRUE)

#Perform PCA, do project onto phyllotactic spiral and retrieve PCA scores

pca <- prcomp(coord[6:35])
pca_scores <- cbind(coord, pca$x)

#Sort datasets so that leaves are arranged top to bottom in the house by their PC1 score within the phyllotactic spiral

sort_PC1 <- pca_scores[order(-pca_scores[36]),]
sort_just_inside <- just_inside[order(-just_inside[4]),]

pca_spiral <- cbind(sort_PC1, sort_just_inside)

#Set a scale (to make leaves porpotionally smaller) and an optional angle delta to rotate leaves

s <- 0.003
delta <- 0

#Create scaled data first

pca_spiral$sx1 <- pca_spiral$x1*s
pca_spiral$sx2 <- pca_spiral$x2*s
pca_spiral$sx3 <- pca_spiral$x3*s
pca_spiral$sx4 <- pca_spiral$x4*s
pca_spiral$sx5 <- pca_spiral$x5*s
pca_spiral$sx6 <- pca_spiral$x6*s
pca_spiral$sx7 <- pca_spiral$x7*s
pca_spiral$sx8 <- pca_spiral$x8*s
pca_spiral$sx9 <- pca_spiral$x9*s
pca_spiral$sx10 <- pca_spiral$x10*s
pca_spiral$sx11 <- pca_spiral$x11*s
pca_spiral$sx12 <- pca_spiral$x12*s
pca_spiral$sx13 <- pca_spiral$x13*s
pca_spiral$sx14 <- pca_spiral$x14*s
pca_spiral$sx15 <- pca_spiral$x15*s

pca_spiral$sy1 <- pca_spiral$y1*s
pca_spiral$sy2 <- pca_spiral$y2*s
pca_spiral$sy3 <- pca_spiral$y3*s
pca_spiral$sy4 <- pca_spiral$y4*s
pca_spiral$sy5 <- pca_spiral$y5*s
pca_spiral$sy6 <- pca_spiral$y6*s
pca_spiral$sy7 <- pca_spiral$y7*s
pca_spiral$sy8 <- pca_spiral$y8*s
pca_spiral$sy9 <- pca_spiral$y9*s
pca_spiral$sy10 <- pca_spiral$y10*s
pca_spiral$sy11 <- pca_spiral$y11*s
pca_spiral$sy12 <- pca_spiral$y12*s
pca_spiral$sy13 <- pca_spiral$y13*s
pca_spiral$sy14 <- pca_spiral$y14*s
pca_spiral$sy15 <- pca_spiral$y15*s

#Using scaled data, create leaves that are rotated and translated to the correct part of the phyllotactic spiral

pca_spiral$nx1 <- (pca_spiral$sx1*cos(pca_spiral$t+delta) - pca_spiral$sy1*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx2 <- (pca_spiral$sx2*cos(pca_spiral$t+delta) - pca_spiral$sy2*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx3 <- (pca_spiral$sx3*cos(pca_spiral$t+delta) - pca_spiral$sy3*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx4 <- (pca_spiral$sx4*cos(pca_spiral$t+delta) - pca_spiral$sy4*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx5 <- (pca_spiral$sx5*cos(pca_spiral$t+delta) - pca_spiral$sy5*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx6 <- (pca_spiral$sx6*cos(pca_spiral$t+delta) - pca_spiral$sy6*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx7 <- (pca_spiral$sx7*cos(pca_spiral$t+delta) - pca_spiral$sy7*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx8 <- (pca_spiral$sx8*cos(pca_spiral$t+delta) - pca_spiral$sy8*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx9 <- (pca_spiral$sx9*cos(pca_spiral$t+delta) - pca_spiral$sy9*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx10 <- (pca_spiral$sx10*cos(pca_spiral$t+delta) - pca_spiral$sy10*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx11 <- (pca_spiral$sx11*cos(pca_spiral$t+delta) - pca_spiral$sy11*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx12 <- (pca_spiral$sx12*cos(pca_spiral$t+delta) - pca_spiral$sy12*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx13 <- (pca_spiral$sx13*cos(pca_spiral$t+delta) - pca_spiral$sy13*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx14 <- (pca_spiral$sx14*cos(pca_spiral$t+delta) - pca_spiral$sy14*sin(pca_spiral$t+delta)) + pca_spiral$x 
pca_spiral$nx15 <- (pca_spiral$sx15*cos(pca_spiral$t+delta) - pca_spiral$sy15*sin(pca_spiral$t+delta)) + pca_spiral$x 

pca_spiral$ny1 <- (pca_spiral$sy1*cos(pca_spiral$t+delta) + pca_spiral$sx1*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny2 <- (pca_spiral$sy2*cos(pca_spiral$t+delta) + pca_spiral$sx2*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny3 <- (pca_spiral$sy3*cos(pca_spiral$t+delta) + pca_spiral$sx3*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny4 <- (pca_spiral$sy4*cos(pca_spiral$t+delta) + pca_spiral$sx4*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny5 <- (pca_spiral$sy5*cos(pca_spiral$t+delta) + pca_spiral$sx5*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny6 <- (pca_spiral$sy6*cos(pca_spiral$t+delta) + pca_spiral$sx6*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny7 <- (pca_spiral$sy7*cos(pca_spiral$t+delta) + pca_spiral$sx7*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny8 <- (pca_spiral$sy8*cos(pca_spiral$t+delta) + pca_spiral$sx8*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny9 <- (pca_spiral$sy9*cos(pca_spiral$t+delta) + pca_spiral$sx9*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny10 <- (pca_spiral$sy10*cos(pca_spiral$t+delta) + pca_spiral$sx10*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny11 <- (pca_spiral$sy11*cos(pca_spiral$t+delta) + pca_spiral$sx11*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny12 <- (pca_spiral$sy12*cos(pca_spiral$t+delta) + pca_spiral$sx12*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny13 <- (pca_spiral$sy13*cos(pca_spiral$t+delta) + pca_spiral$sx13*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny14 <- (pca_spiral$sy14*cos(pca_spiral$t+delta) + pca_spiral$sx14*sin(pca_spiral$t+delta)) + pca_spiral$y 
pca_spiral$ny15 <- (pca_spiral$sy15*cos(pca_spiral$t+delta) + pca_spiral$sx15*sin(pca_spiral$t+delta)) + pca_spiral$y 

#Visualize!!! Using geom_segment()

size=0.1

p <- ggplot()
p + geom_segment(data=pca_spiral, aes(x=nx1, y=ny1, xend=nx7, yend=ny7), size=size) +
geom_segment(data=pca_spiral, aes(x=nx7, y=ny7, xend=nx8, yend=ny8), size=size) +
geom_segment(data=pca_spiral, aes(x=nx8, y=ny8, xend=nx9, yend=ny9), size=size) +
geom_segment(data=pca_spiral, aes(x=nx9, y=ny9, xend=nx10, yend=ny10), size=size) +
geom_segment(data=pca_spiral, aes(x=nx10, y=ny10, xend=nx11, yend=ny11), size=size) +
geom_segment(data=pca_spiral, aes(x=nx11, y=ny11, xend=nx12, yend=ny12), size=size) +
geom_segment(data=pca_spiral, aes(x=nx12, y=ny12, xend=nx13, yend=ny13), size=size) +
geom_segment(data=pca_spiral, aes(x=nx13, y=ny13, xend=nx14, yend=ny14), size=size) +
geom_segment(data=pca_spiral, aes(x=nx14, y=ny14, xend=nx15, yend=ny15), size=size) + 
geom_segment(data=pca_spiral, aes(x=nx15, y=ny15, xend=nx6, yend=ny6), size=size) + 
geom_segment(data=pca_spiral, aes(x=nx6, y=ny6, xend=nx1, yend=ny1), size=size) + 
geom_polygon(data=data, aes(x=new_x, y=new_y), colour="black", fill="#00000000") + coord_fixed() + 
theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())


