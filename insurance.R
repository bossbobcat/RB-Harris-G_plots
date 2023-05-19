## Add libraries
library(tidyverse)
library(ggfortify)
library(survival)
library(patchwork)
library(ggsurvfit)



##########
## Data ##
##########

x=sort(c(08.0, 18.2, 22.7, 16.3, 11.1, 26.9, 26.6, 18.9, 24.3, 20.6,
         40.7, 56.9, 49.3, 55.0, 49.4, 50.7, 46.6, 87.9, 79.2, 78.9,
         82.1, 63.7, 61.7, 68.8, 60.8, 69.2, 55.9, 54.8, 65.4, 53.7,
         65.5, 52.7, 52.9, 50.7, 61.3, 49.6, 47.5, 58.9, 47.4, 56.0,
         46.9, 46.5, 57.9, 45.7, 44.5, 53.1, 44.1, 41.8, 48.2, 37.1,
         32.7, 37.6, 42.8, 47.4, 35.6, 32.2, 30.1, 31.2))

## Sorted Data
df = data.frame(x=x,evt=1)
km = survfit(Surv(x,evt)~1,data=df)

## Survival data
y_surv= km$surv



############################################
## PLOT 1: Survival Function of the Model ##
############################################



rbhw <-function(x,delta,v,theta,lambda)
{
  pgamma(-log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)), delta)
}

## survival model #############
y_main=rbhw(x,
            delta = 9.8921e-01, 
            v = 3.5825e-01,
            theta =  6.2749e+02, 
            lambda =6.9285e-01)

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
##ggsave("survival_insurance.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



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
##ggsave("ecdf_insurance.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



##################################
## PLOT 3: Scaled TTT Transform ##
##################################



df = tibble(x=(1:length(xs))/length(xs) , f = (cumsum(xs) + (length(xs)-(1:length(x)))*xs)/ sum(xs))
gg = df %>%
  ggplot(aes(x=x,y=f))+
  geom_line(linewidth=1)+
  geom_point(col="green",size=2)+
  geom_abline(intercept=0.1704504, slope= 0.8295496,linetype=2)+
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
##ggsave("TTT_insurance.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



##############################
## PLOT 4: HRF of the model ##
##############################



x=sort(x)
hrf.rbhw = function(x,delta,v,theta,lambda){
  ((1/gamma(delta))*(( -log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)))^(delta-1))*((theta^(1/v))*lambda*(x^(lambda-1))*exp(-x^lambda))/((1-(1-theta)*(exp(-x^lambda))^v)^(1+(1/v)))
   
  )/(pgamma(-log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)), delta))
}

hrf.df=tibble(x=x , y = hrf.rbhw(x,
                                 delta =9.8921e-01, 
                                 v = 3.5825e-01,
                                 theta =  6.2749e+02, 
                                 lambda =6.9285e-01))

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
##ggsave("h_insurance.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



########################################
## PLOT 5: Expected Probability Plots ##
########################################



F_observed = ((1:length(x))-0.375)/(length(x)+0.25)
F.RBHW = 1- rbhw(x,
                 delta =9.8921e-01, 
                 v = 3.5825e-01,
                 theta =  6.2749e+02, 
                 lambda =6.9285e-01)

SS.RBHW = sum((F.RBHW - F_observed)^2)



#GELLoG

GELLoG_cdf = function(x,lambda,c,alpha, delta) {
  (pgamma(-log(1-((1-((1+lambda+lambda*x)/(1+lambda))*(exp(-lambda*x)/(1+x^c)))^alpha)),delta))
}

F.GELLoG = GELLoG_cdf(x,lambda=0.2317,c=0.00000052945,alpha=0.0039148, delta=14.8790)

SS.GELLoG = sum((F.GELLoG - F_observed)^2)    

### APTLW 

APTLW_cdf <- function(x,theta,alpha, beta ,lambda) {
  ((alpha^((1-exp(-2*lambda*x^beta))^theta))-1)/(alpha-1)
}

F.APTLW = APTLW_cdf(x,theta=  0.5283, alpha=  104.97, beta=2.1481  ,lambda=0.000174117)
SS.APTLW = sum((F.APTLW - F_observed)^2)   

#TLGW

TLGW_cdf <- function(data,alpha, theta, lambda,beta) {
  ((1-exp(-(lambda*x)^beta))^(theta*alpha))*(2-(1-exp(-(lambda*x)^beta))^theta)^(alpha) 
}

F.TLGW = TLGW_cdf(x,alpha=94.941,
                  theta=  0.013352 ,
                  lambda= 0.010763,
                  beta=9.1307)
SS.TLGW = sum((F.TLGW - F_observed)^2)   


#KW
KW_cdf <- function(x,a,b,c, lambda){
  1-(1-(1-exp(-(lambda*x)^c))^a)^b
}

F.KW = KW_cdf(x,
              a = 192.45,
              b = 2814.5,
              lambda = 1068.1,
              c = 0.1067)

SS.KW = sum((F.KW - F_observed)^2)

## GLLoGW
GLLoGW_cdf <- function(x,c,beta,delta,theta){
  1-pgamma(-(theta^(-1))*log(1-((1+x^c)^(-1))*exp(-1*x^beta)), delta)
}

F.GLLoGW = GLLoGW_cdf(x,
                      c = 0.8350, 
                      beta = 0.1109,
                      delta = 4.2111, 
                      theta = 0.0024)

SS.GLLoGW = sum((F.GLLoGW - F_observed)^2)

## group data
df.cdf = tibble(F_observed,F.RBHW,F.APTLW,F.GELLoG,F.TLGW,F.KW,F.GLLoGW) %>%
  pivot_longer(cols=-F_observed, names_to = "dist", values_to = "y")

gg= df.cdf %>%
  ggplot(aes(x=F_observed,y=y, col=dist))+
  geom_line()+
  scale_colour_manual(name="", #legend name with a title "Parameters" if you want
                      labels = c(paste0("(SS=",round(SS.APTLW,4),") APTLW"),
                                 paste0("(SS=",round(SS.GELLoG,4),") GELLoG"),
                                 paste0("(SS=",round(SS.GELLoG,4),") GLLoGW"),
                                 paste0("(SS=",round(SS.KW,4),") KW"),
                                 paste0("(SS=",round(SS.RBHW,4),") RB-Harris-W"),
                                 paste0("(SS=",round(SS.TLGW,4),") TLGW")
                      ),
                      values = c("red", "magenta","orange","green","black","blue"))+
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
##ggsave("pp_insurance.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution