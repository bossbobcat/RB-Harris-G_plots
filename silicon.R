## Add libraries
library(tidyverse)
library(ggfortify)
library(survival)
library(patchwork)
library(ggsurvfit)

##########################
## Silicon Nitride Data ##
##########################

x=sort(c(5.50, 5.00, 4.90, 6.40, 5.10, 5.20, 5.20, 5.00, 4.70, 4.00,
         4.50, 4.20, 4.10, 4.56, 5.01, 4.70, 3.13, 3.12, 2.68, 2.77,
         2.70, 2.36, 4.38, 5.73, 4.35, 6.81, 1.91, 2.66, 2.61, 1.68,
         2.04, 2.08, 2.13, 3.80, 3.73, 3.71, 3.28, 3.90, 4.00, 3.80,
         4.10, 3.90, 4.05, 4.00, 3.95, 4.00, 4.50, 4.50, 4.20, 4.55,
         4.65, 4.10, 4.25, 4.30, 4.50, 4.70, 5.15, 4.30, 4.50, 4.90,
         5.00, 5.35, 5.15, 5.25, 5.80, 5.85, 5.90, 5.75, 6.25, 6.05,
         5.90, 3.60, 4.10, 4.50, 5.30, 4.85, 5.30, 5.45, 5.10, 5.30,
         5.20, 5.30, 5.25, 4.75, 4.50, 4.20, 4.00, 4.15, 4.25, 4.30,
         3.75, 3.95, 3.51, 4.13, 5.40, 5.00, 2.10, 4.60, 3.20, 2.50,
         4.10, 3.50, 3.20, 3.30, 4.60, 4.30, 4.30, 4.50, 5.50, 4.60,
         4.90, 4.30, 3.00, 3.40, 3.70, 4.40, 4.90, 4.90, 5.00))

df = data.frame(x=x,evt=1)

km = survfit(Surv(x,evt)~1,data=df)

## Survival data
y_surv= km$surv



############################################
## PLOT 1: Survival Function of the Model ##
############################################



rbhw <-function(x,delta,v,theta,lambda){
  pgamma(-log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)), delta)
}

## survival model #############
y_main=rbhw(x,
            delta = 0.4927957, 
            v =  0.2844730, 
            theta = 57.4990599, 
            lambda = 1.6975765)

## group everything
df$y_main = y_main
df.surv = data.frame(x=unique(x),y_surv=y_surv)

df.grp = left_join(df,df.surv)

gg = df.grp %>%
  ggplot(aes(x=x))+
  geom_line(aes(y=y_main,col=I("black")),size=1) +
  geom_step(aes(y=y_surv,col="red"),size=1)+
  scale_colour_manual(name=NULL, #legend name
                      values = c("black","red"), #specify and change colors of all lines
                      labels = c("RB-Harris-W","Kaplan-Meier estimate"))+
  
  labs(y="Survival Function", x = "x")+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 15,color="black"), #axis font 
        axis.title = element_text(size=15,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        #plot.margin = margin(1, 1, 1, 1, "cm"), # t, r, b, l
        legend.position = c(1,1), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1,1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
##ggsave("survival_silicon.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



##################
## PLOT 2: ECDF ##
##################



gg = df.grp %>%
  ggplot(aes(x=x))+
  geom_line(aes(y=1-y_main,col=I("black")),size=1) +
  geom_step(aes(y=1-y_surv,col="red"),size=1)+
  scale_colour_manual(name=NULL, #legend name
                      values = c("black","red"), #specify and change colors of all lines
                      labels = c("RB-Harris-W","ECDF"))+
  
  labs(y="F(x)", x = "x")+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 15,color="black"), #axis font 
        axis.title = element_text(size=15,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        #plot.margin = margin(1, 1, 1, 1, "cm"), # t, r, b, l
        legend.position = c(0.95,0.2), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1,1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
##ggsave("ecdf_silicon.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



##################################
## PLOT 3: Scaled TTT Transform ##
##################################



df = tibble(x=(1:length(xs))/length(xs) , f = (cumsum(xs) + (length(xs)-(1:length(x)))*xs)/ sum(xs))
gg = df %>%
  ggplot(aes(x=x,y=f))+
  geom_line(linewidth=1)+
  geom_point(col="green",size=2)+
  geom_abline(intercept=0.3884053, slope=0.6115947,linetype=2)+
  labs(y="Scaled TTT-Transform")+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 15,color="black"), #axis font 
        axis.title = element_text(size=15,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        #plot.margin = margin(1, 1, 1, 1, "cm"), # t, r, b, l
        legend.position = c(1,1), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1,1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
##ggsave("TTT_silicon.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



##############################
## PLOT 4: HRF of the model ##
##############################



hrf.rbhw = function(x,delta,v,theta,lambda){
  ((1/gamma(delta))*(( -log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)))^(delta-1))*((theta^(1/v))*lambda*(x^(lambda-1))*exp(-x^lambda))/((1-(1-theta)*(exp(-x^lambda))^v)^(1+(1/v)))
   
  )/(pgamma(-log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)), delta))
}

hrf.df = tibble(x=x, y = hrf.rbhw(x,
                                  delta = 0.4927957, 
                                  v =  0.2844730, 
                                  theta = 57.4990599, 
                                  lambda = 1.6975765))

gg= hrf.df %>%
  ggplot(aes(x=x,y=y))+
  geom_line(linewidth=1,col="black")+
  labs(y="h(x) of RB-Harris-W")+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 15,color="black"), #axis font 
        axis.title = element_text(size=15,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        #plot.margin = margin(1, 1, 1, 1, "cm"), # t, r, b, l
        legend.position = c(1,1), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1,1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
##ggsave("h_silicon.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



###############################
## PLOT 5: Probability Plots ##
###############################



## RB-Harris-W
F_observed = ((1:length(x))-0.375)/(length(x)+0.25)

F.RBHW = 1- rbhw(x,
                 delta = 0.4927957, 
                 v = 0.2844730, 
                 theta = 57.4990599, 
                 lambda = 1.6975765)

SS.RBHW = sum((F.RBHW - F_observed)^2)

## GELLoG
GELLoG_cdf = function(x,lambda,c,alpha, delta) {
  (pgamma(-log(1-((1-((1+lambda+lambda*x)/(1+lambda))*(exp(-lambda*x)/(1+x^c)))^alpha)),delta))
}

F.GELLoG = GELLoG_cdf(x,
                      lambda = 3.7014,
                      c = 1.0472,
                      alpha = 0.4907, 
                      delta = 15.8740)

SS.GELLoG = sum((F.GELLoG - F_observed)^2)    

## APTLW 
APTLW_cdf <- function(x,theta,alpha,beta,lambda) {
  ((alpha^((1-exp(-2*lambda*x^beta))^theta))-1)/(alpha-1)
}

F.APTLW = APTLW_cdf(x,
                    theta = 0.7681,
                    alpha = 60.4077, 
                    beta = 3.3023,
                    lambda = 0.0059)

SS.APTLW = sum((F.APTLW - F_observed)^2)   

## TLGW
TLGW_cdf <- function(data,alpha, theta, lambda,beta) {
  ((1-exp(-(lambda*x)^beta))^(theta*alpha))*(2-(1-exp(-(lambda*x)^beta))^theta)^(alpha) 
}

F.TLGW = TLGW_cdf(x,
                  alpha = 101.02,
                  theta =  0.012318 ,
                  lambda = 0.15012,
                  beta = 16.077)

SS.TLGW = sum((F.TLGW - F_observed)^2)   

## KW
KW_cdf <- function(x,a,b,c, lambda){
  1-(1-(1-exp(-(lambda*x)^c))^a)^b
}

F.KW = KW_cdf(x,
              a=0.8261,
              b=0.5374,
              lambda=0.2326,
              c=5.4234)

SS.KW = sum((F.KW - F_observed)^2)

## TODO: Add GLLoGW here

## group data
df.cdf = tibble(F_observed,F.RBHW,F.APTLW,F.GELLoG,F.TLGW,F.KW) %>%
  pivot_longer(cols=-F_observed, names_to = "dist", values_to = "y")

gg= df.cdf %>%
  ggplot(aes(x=F_observed,y=y, col=dist))+
  geom_line()+
  scale_colour_manual(name="", #legend name with a title "Parameters" if you want
                      labels = c(paste0("(SS=",round(SS.APTLW,4),") APTLW"),
                                 paste0("(SS=",round(SS.GELLoG,4),") GELLoG"),
                                 paste0("(SS=",round(SS.KW,4),") KW"),
                                 paste0("(SS=",round(SS.RBHW,4),") RB-Harris-W"),
                                 paste0("(SS=",round(SS.TLGW,4),") TLGW")
                      ),
                      values = c("red", "magenta","green","black","blue"))+
  geom_line(aes(x=F_observed,y=F_observed))+
  theme_bw()+
  labs(x="Observed Probability", y="Expected Probability")+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font 
        axis.title = element_text(size=20,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 18),
        legend.position = c(0.95,0.3), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
##ggsave("pp_silicon.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution