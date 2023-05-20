library(tidyverse)
######### pdfs B-Harris-SHL ###########
pdfs=function(x,delta,v,theta){
  y=((1/gamma(delta))*(theta^(1/v))*((2*exp(-x))/((1+exp(-x))^2))*((-log(1-((theta*(1-(1-exp(-x))/(1+exp(-x)))^v)/(1-(1-theta)*(1-(1-exp(-x))/(1+exp(-x)))^v))^(1/v)))^(delta-1))/(1-(1-theta)*(1-(1-exp(-x))/(1+exp(-x)))^v)^(1+(1/v))
  )
  return(y)
}

x_sim=seq(0.00001,1.5,by=0.001)
densities = data.frame(x=x_sim,
                 y1=pdfs(x_sim,0.3,0.09,.04),
                 y2=pdfs(x_sim,2.9,.3,.6),
                 y3=pdfs(x_sim,1.9,9.9,69.9),
                 y4=pdfs(x_sim,4.8,12.8,90.0),
                 y5=pdfs(x_sim,3.0,10.0,78.0)) %>%
pivot_longer(cols = -x,names_to = "Densities", values_to = "y")
  


densities %>% 
  ggplot(aes(x=x,y=y,colour=Densities)) +
  geom_line(lwd=.8) +
  scale_colour_manual("",
                      values = c("red", "magenta", "green","cyan","black","blue"),
                      labels = c(expression(delta ~ "=0.3, v=0.09,"~ theta ~"=0.04"),
                                 expression(delta ~ "=2.9, v=0.3,"~ theta ~"=0.6"),
                                 expression(delta ~ "=1.9, v=9.9,"~ theta ~"=69.9"),
                                 expression(delta ~ "=4.8, v=12.8,"~ theta ~"=90.0"),
                                 expression(delta ~ "=3.0, v=10.0,"~ theta ~"=78.80")))+
  scale_linetype_manual("",
                        values = c(1:5),
                        labels = c(expression(delta ~ "=0.3, v=0.09,"~ theta ~"=0.04"),
                                   expression(delta ~ "=2.9, v=0.3,"~ theta ~"=0.6"),
                                   expression(delta ~ "=1.9, v=9.9,"~ theta ~"=69.9"),
                                   expression(delta ~ "=4.8, v=12.8,"~ theta ~"=90.0"),
                                   expression(delta ~ "=3.0, v=10.0,"~ theta ~"=78.80"))) +
  ylim(0,4)+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 15,color="black"), #axis font 
        axis.title = element_text(size=15,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        #legend.title = element_text(face="bold",size = 15),
        legend.text.align = 0,
        legend.position = c(0.99,1.0), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        #legend.spacing.x = unit(1.0, 'cm'),
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) + #make transparent background in legend
  labs(y=bquote(f[RB-Harris-SHL](x))) -> gg1
print(gg1)

ggsave("pdfs_SHL.eps", gg1, width=18, height=17, units="cm", dpi=1080) #save with a specific dimension and resolution







######### hrts B-Harris-SHL ###########
hrt=function(x,delta,v,theta){
  y=((1/gamma(delta))*(theta^(1/v))*((2*exp(-x))/((1+exp(-x))^2))*((-log(1-((theta*(1-(1-exp(-x))/(1+exp(-x)))^v)/(1-(1-theta)*(1-(1-exp(-x))/(1+exp(-x)))^v))^(1/v)))^(delta-1))/(1-(1-theta)*(1-(1-exp(-x))/(1+exp(-x)))^v)^(1+(1/v))
  )/(pgamma(-log(1-((theta*(1-(1-exp(-x))/(1+exp(-x)))^v)/(1-(1-theta)*(1-(1-exp(-x))/(1+exp(-x)))^v))^(1/v)),delta)
  )
  return(y)
}

x_sim=seq(0.00001,1.5,by=0.001)
hratios = data.frame(x=x_sim,
                     y1=hrt(x_sim,0.2,0.1,.04),
                     y2=hrt(x_sim,1.9,4.9,34.2),
                     y3=hrt(x_sim,3.9,3.0,7.2),
                     y4=hrt(x_sim,2.8,.01,5.0),
                     y5=hrt(x_sim,2.6,29.0,3.9)) %>%
  pivot_longer(cols = -x,names_to = "Hratios", values_to = "y")



hratios %>% 
  ggplot(aes(x=x,y=y,colour=Hratios)) +
  geom_line(lwd=.8) +
  scale_colour_manual("",
                      values = c("red", "magenta", "green","cyan","black","blue"),
                      labels = c(expression(delta ~ "=0.2, v=0.1,"~ theta ~"=0.04"),
                                 expression(delta ~ "=1.9, v=4.9,"~ theta ~"=34.2"),
                                 expression(delta ~ "=3.9, v=3.0,"~ theta ~"=7.2"),
                                 expression(delta ~ "=2.8, v=0.01,"~ theta ~"=5.0"),
                                 expression(delta ~ "=2.6, v=29.0,"~ theta ~"=3.9")))+
  scale_linetype_manual("",
                        values = c(1:5),
                        labels = c(expression(delta ~ "=0.2, v=0.1,"~ theta ~"=0.04"),
                                   expression(delta ~ "=1.9, v=4.9,"~ theta ~"=34.2"),
                                   expression(delta ~ "=3.9, v=3.0,"~ theta ~"=7.2"),
                                   expression(delta ~ "=2.8, v=0.01,"~ theta ~"=5.0"),
                                   expression(delta ~ "=2.6, v=29.0,"~ theta ~"=3.9"))) +
  ylim(0,5)+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 15,color="black"), #axis font 
        axis.title = element_text(size=15,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        #legend.title = element_text(face="bold",size = 15),
        legend.text.align = 0,
        legend.position = c(0.99,1.0), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        #legend.spacing.x = unit(1.0, 'cm'),
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) + #make transparent background in legend
  labs(y=bquote(h[RB-Harris-SHL](x))) -> gg2
print(gg2)

ggsave("hrfs_SHL.eps", gg2, width=18, height=17, units="cm", dpi=1080) #save with a specific dimension and resolution


library(ggpubr)
gg= ggarrange(gg1,gg2,
              ncol = 1,
          labels = "AUTO"
            )
print(gg)
ggsave("pdfhrf_SHL.eps", gg, width=18, height=25, units="cm", dpi=1080) #save with a specific dimension and resolution

