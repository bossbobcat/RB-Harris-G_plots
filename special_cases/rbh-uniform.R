library(tidyverse)
######### pdfs B-Harris-Weibull ###########
f1=function(x,delta,v,theta){
  y=((1/gamma(delta))*(( -log(1-((theta*((x-min(x))/(max(x)-min(x)))^v)/(1-(1-theta)*((x-min(x))/(max(x)-min(x)))^v))^(1/v)))^(delta-1))*((theta^(1/v))*(1/(max(x)-min(x)))/((1-(1-theta)*((x-min(x))/(max(x)-min(x))))^v)^(1+(1/v)))
  )
  return(y)
}

f2=function(x){
  y=dunif(x,min(x),max(x))
  return(y)
}

x_sim=seq(0,1.5,by=0.001)
densities = data.frame(x=x_sim,
                       y1=f2(x_sim),
                       y1 =f1(x_sim,2.0,0.9,1.0,.1),
                       y2=f1(x_sim,0.7,0.8,0.9),
                       y3=f1(x_sim,1.9,2.4,0.4),
                       y4=f1(x_sim,0.6,0.3,0.2),
                       y5=f1(x_sim,0.1,1.9,0.5)
) %>%
pivot_longer(cols = -x,names_to = "Densities", values_to = "y")



densities %>%
  ggplot(aes(x=x,y=y,colour=Densities)) +
  geom_line(lwd=0.8) +
  scale_colour_manual("",
                      values = c("red", "magenta", "green","cyan","black","blue"),
                      labels = c(expression(delta ~ "=2.0, v=0.9,"~ theta ~"=1.0," ~lambda~ "=0.1"),
                                 expression(delta ~ "=0.7, v=0.8,"~ theta ~"=0.9," ~lambda~ "=3.8"),
                                 expression(delta ~ "=1.9, v=2.4,"~ theta ~"=0.4," ~lambda~ "=2.2"),
                                 expression(delta ~ "=0.6, v=0.3,"~ theta ~"=0.2," ~lambda~ "=3.9"),
                                 expression(delta ~ "=0.1, v=1.9,"~ theta ~"=0.5," ~lambda~ "=4.8")))+
  scale_linetype_manual("",
                        values = c(1:5),
                        labels = c(expression(delta ~ "=2.0, v=0.9,"~ theta ~"=1.0," ~lambda~ "=0.1"),
                                   expression(delta ~ "=0.7, v=0.8,"~ theta ~"=0.9," ~lambda~ "=3.8"),
                                   expression(delta ~ "=1.9, v=2.4,"~ theta ~"=0.4," ~lambda~ "=2.2"),
                                   expression(delta ~ "=0.6, v=0.3,"~ theta ~"=0.2," ~lambda~ "=3.9"),
                                   expression(delta ~ "=0.1, v=1.9,"~ theta ~"=0.5," ~lambda~ "=4.8"))) +
  ylim(0,3)+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font
        axis.title = element_text(size=20,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 16),  #legend font
        #legend.title = element_text(face="bold",size = 15),
        legend.position = c(0.99,1), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        legend.spacing.x = unit(1.0, 'cm'),
        axis.line.x.bottom = element_line(color = "black",size = 0.7), #control axis and its color
        axis.line.y.left = element_line(color = "black", size = 0.7), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) + #make transparent background in legend
  labs(y=bquote(f[RB-Harris-W](x))) -> gg1
print(gg1)

#ggsave("pdfs_w.eps", gg1, width=18, height=17, units="cm", dpi=1080) #save with a specific dimension and resolution







######### hrts B-Harris-Weibull ###########
hrt=function(x,delta,v,theta,lambda){
  y=((1/gamma(delta))*(( -log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)))^(delta-1))*((theta^(1/v))*lambda*(x^(lambda-1))*exp(-x^lambda))/((1-(1-theta)*(exp(-x^lambda))^v)^(1+(1/v)))
  )/(pgamma(-log(1-((theta*(exp(-x^lambda))^v)/(1-(1-theta)*(exp(-x^lambda))^v))^(1/v)), delta)
  )
  return(y)
}

x_sim=seq(0,1.5,by=0.001)
hratios = data.frame(x=x_sim,
                     y1=hrt(x_sim,2.0,0.9,1.0,0.1),
                     y2=hrt(x_sim,0.7,0.8,0.9,3.8),
                     y3=hrt(x_sim,1.9,2.4,0.4,2.2),
                     y4=hrt(x_sim,0.6,0.3,0.2,3.9),
                     y5=hrt(x_sim,0.1,1.9,0.5,4.8)) %>%
  pivot_longer(cols = -x,names_to = "Hratios", values_to = "y")


hratios %>%
  ggplot(aes(x=x,y=y,colour=Hratios)) +
  geom_line(lwd=0.8) +
  scale_colour_manual("",
                      values = c("red", "magenta", "green","cyan","black","blue"),
                      labels =  c(expression(delta ~ "=2.0, v=0.9,"~ theta ~"=1.0," ~lambda~ "=0.1"),
                                  expression(delta ~ "=0.7, v=0.8,"~ theta ~"=0.9," ~lambda~ "=3.8"),
                                  expression(delta ~ "=1.9, v=2.4,"~ theta ~"=0.4," ~lambda~ "=2.2"),
                                  expression(delta ~ "=0.6, v=0.3,"~ theta ~"=0.2," ~lambda~ "=3.9"),
                                  expression(delta ~ "=0.1, v=1.9,"~ theta ~"=0.5," ~lambda~ "=4.8")))+
  scale_linetype_manual("",
                        values = c(1:5),
                        labels =  c(expression(delta ~ "=2.0, v=0.9,"~ theta ~"=1.0," ~lambda~ "=0.1"),
                                    expression(delta ~ "=0.7, v=0.8,"~ theta ~"=0.9," ~lambda~ "=3.8"),
                                    expression(delta ~ "=1.9, v=2.4,"~ theta ~"=0.4," ~lambda~ "=2.2"),
                                    expression(delta ~ "=0.6, v=0.3,"~ theta ~"=0.2," ~lambda~ "=3.9"),
                                    expression(delta ~ "=0.1, v=1.9,"~ theta ~"=0.5," ~lambda~ "=4.8"))) +
  ylim(0,12)+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font
        axis.title = element_text(size=20,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 16),  #legend font
        #legend.title = element_text(face="bold",size = 15),
        legend.position = c(0.99,1), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        legend.spacing.x = unit(1.0, 'cm'),
        axis.line.x.bottom = element_line(color = "black",size = 0.7), #control axis and its color
        axis.line.y.left = element_line(color = "black", size = 0.7), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) + #make transparent background in legend
  labs(y=bquote(h[RB-Harris-W](x))) -> gg2
print(gg2)

#ggsave("hrfs_RBH_W.eps", gg2, width=18, height=17, units="cm", dpi=1080) #save with a specific dimension and resolution


library(ggpubr)
gg= ggarrange(gg1,gg2,
              ncol = 1,
          labels = "AUTO"
            )
print(gg)
#ggsave("pdfhrf_Weibull.eps", gg, width=15, height=25, units="cm", dpi=1080) #save with a specific dimension and resolution

