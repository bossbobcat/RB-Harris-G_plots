### Code plotting

## Add libraries
library(tidyverse)
library(ggfortify)
library(survival)
library(patchwork)
#library(survminer)
library(ggsurvfit)


## Data 

x=c(4.4130, 3.0525, 4.6955, 7.4810, 5.1915, 3.6335, 6.6100, 8.2490,
    5.8325, 3.0075, 5.4275, 3.0610, 3.3280, 1.7200, 2.9270, 5.3425, 5.0175, 2.6210, 2.1720, 2.5715,
    3.8150, 7.3020, 3.9515, 3.1850, 1.7685, 3.1635, 2.3650, 1.6075, 4.6420, 6.4390, 4.4065, 5.0215,
    3.6300, 2.9925, 3.2060, 1.6975, 2.2120, 4.9675, 3.9200, 4.7750, 1.7495, 1.8755, 3.4840, 1.6430,
    5.0790, 4.0540, 3.3485, 3.5755, 3.2800, 1.0385, 1.8890, 1.4940, 1.6680, 3.4070, 4.1625, 3.9270,
    4.2755, 1.6140, 3.7430, 3.3125, 3.0700, 2.4545, 2.3305, 2.6960, 6.0210, 4.3480, 0.9075, 1.6635,
    2.7030, 3.0910, 0.5205, 0.9000, 2.4745, 2.0445, 1.6795, 1.0350, 1.6490, 2.6585, 2.7210, 2.2785,
    2.1460, 1.2500, 3.2675, 2.3240, 2.3485, 2.7295, 2.0600, 1.9610, 1.6095, 0.7010, 1.2190, 1.6285,
    1.8160, 1.6165, 1.5135, 1.1760, 0.6025, 1.6090, 1.4630, 1.3005, 1.0325, 1.5145, 1.0290, 1.1630,
    1.2530, 0.9615)

## Sorted Data
x = sort(x)
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
            delta =0.191359,   
            v =13.524331, 
            theta =115.349737,  
            lambda =1.555517)

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
##ggsave("survival_covid.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



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
##ggsave("ecdf_covid.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



##################################
## PLOT 3: Scaled TTT Transform ##
##################################



xs = sort(x)
df = tibble(x=(1:length(xs))/length(xs) , f = (cumsum(xs) + (length(xs)-(1:length(x)))*xs)/ sum(xs))
gg = df %>%
  ggplot(aes(x=x,y=f))+
  geom_line(linewidth=1)+
  geom_point(col="green",size=2)+
  geom_abline(intercept=0.179, slope=0.821,linetype=2)+
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
##ggsave("TTT_covid.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



##############################
## PLOT 4: HRF of the model ##
##############################



x=sort(x)
hrf.rbhw = function(x,delta,v,theta,lambda){
  ((1/gamma(delta))*(( -log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)))^(delta-1))*((theta^(1/v))*lambda*(x^(lambda-1))*exp(-x^lambda))/((1-(1-theta)*(exp(-x^lambda))^v)^(1+(1/v)))
   
  )/(pgamma(-log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)), delta))
}

hrf.df=tibble(x=x , y = hrf.rbhw(x,
                                 delta =0.191359,   
                                 v =13.524331, 
                                 theta =115.349737,  
                                 lambda =1.555517))

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
##ggsave("h_covid.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



########################################
## PLOT 5: Expected Probability Plots ##
########################################



F_observed = ((1:length(x))-0.375)/(length(x)+0.25)
F.RBHW = 1- rbhw(x,
                 delta =0.191359,   
                 v =13.524331, 
                 theta =115.349737,  
                 lambda =1.555517)

SS.RBHW = sum((F.RBHW - F_observed)^2)



#GELLoG

GELLoG_cdf = function(x,lambda,c,alpha, delta) {
  (pgamma(-log(1-((1-((1+lambda+lambda*x)/(1+lambda))*(exp(-lambda*x)/(1+x^c)))^alpha)),delta))
}

F.GELLoG = GELLoG_cdf(x,lambda=0.9081,c=1.3517,alpha=1.4120, delta=3.0399)

SS.GELLoG = sum((F.GELLoG - F_observed)^2)    

### APTLW 

APTLW_cdf <- function(x,theta,alpha, beta ,lambda) {
  ((alpha^((1-exp(-2*lambda*x^beta))^theta))-1)/(alpha-1)
}

F.APTLW = APTLW_cdf(x,theta=  5.8553, alpha=0.000027609  , beta=0.5942  ,lambda=0.2816)
SS.APTLW = sum((F.APTLW - F_observed)^2)   

#TLGW

TLGW_cdf <- function(data,alpha, theta, lambda,beta) {
  ((1-exp(-(lambda*x)^beta))^(theta*alpha))*(2-(1-exp(-(lambda*x)^beta))^theta)^(alpha) 
}

F.TLGW = TLGW_cdf(x,alpha=2.2614,
                  theta=1.6008 ,
                  lambda= 0.3910,
                  beta=0.9870)
SS.TLGW = sum((F.TLGW - F_observed)^2)   


#KW
KW_cdf <- function(x,a,b,c, lambda){
  1-(1-(1-exp(-(lambda*x)^c))^a)^b
}

F.KW = KW_cdf(x,
              a = 23.7652,
              b = 526.5723,
              ## Values for c and lambda are flipped in the paper.
              ## This configuration results in smaller SS.KW.
              c = 0.1835,  
              lambda = 2.4539)

SS.KW = sum((F.KW - F_observed)^2)

## GLLoGW
GLLoGW_cdf <- function(x,c,beta,delta,theta){
  1-pgamma(-(theta^(-1))*log(1-((1+x^c)^(-1))*exp(-1*x^beta)), delta)
}

F.GLLoGW = GLLoGW_cdf(x,
                      c = 11.9332, 
                      beta = 0.9698,
                      delta = 0.0786, 
                      theta = 0.4962)

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
##ggsave("pp_covid.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution
