library(tidyverse)
######### pdfs B-Harris-Log-Logistic ###########
f1=function(x,delta,v,theta,c){
  y=((1/gamma(delta))*(theta^(1/v))*(c*(x^(c-1))*(1+x^c)^(-2))*((-log(1-((theta*(1+x^c)^(-v))/(1-(1-theta)*(1+x^c)^(-v)))^(1/v)))^(delta-1))/((1-(1-theta)*(1+x^c)^(-v)))^(1+(1/v))
     
  )
  return(y)
}

x_sim=seq(0.01,1.5,by=0.001)
densities = data.frame(x=x_sim,
                 y1=f1(x_sim,0.6,0.4,0.3,0.5),
                 y2=f1(x_sim,3.0,0.7,1.0,2.5),
                 y3=f1(x_sim,1.0,0.2,0.9,6.8),
                 y4=f1(x_sim,1.8,1.9,.8,5.4),
                 y5=f1(x_sim,.9,1.8,19.8,6.8)) %>%
pivot_longer(cols = -x,names_to = "Densities", values_to = "y")
  


densities %>% 
  ggplot(aes(x=x,y=y,colour=Densities)) +
  geom_line(lwd=1) +
  scale_colour_manual("",
                      values = c("red", "magenta", "green","cyan","black","blue"),
                      labels = c(expression(delta ~ "=0.6, v =0.4,"~ theta ~"=0.3, c =0.5"),
                                 expression(delta ~ "=3.0, v =0.7,"~ theta ~"=1.0, c =2.5"),
                                 expression(delta ~ "=1.0, v =0.2,"~ theta ~"=0.9, c =6.8"),
                                 expression(delta ~ "=1.8, v =1.9,"~ theta ~"=0.8, c =5.4"),
                                 expression(delta ~ "=0.9, v =1.8,"~ theta ~"=19.8, c =6.8")))+
  scale_linetype_manual("",
                        values = c(1:5),
                        labels = c(expression(delta ~ "=0.6, v =0.4,"~ theta ~"=0.3, c =0.5"),
                                   expression(delta ~ "=3.0, v =0.7,"~ theta ~"=1.0, c =2.5"),
                                   expression(delta ~ "=1.0, v =0.2,"~ theta ~"=0.9, c =6.8"),
                                   expression(delta ~ "=1.8, v =1.9,"~ theta ~"=0.8, c =5.4"),
                                   expression(delta ~ "=0.9, v =1.8,"~ theta ~"=19.8, c =6.8"))) +
  ylim(0,4)+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 15,color="black"), #axis font 
        axis.title = element_text(size=15,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        #legend.title = element_text(face="bold",size = 15),
        legend.text.align = 0,
        legend.position = c(0.975,1.0), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "left",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        #legend.spacing.x = unit(1.0, 'cm'),
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) + #make transparent background in legend
  labs(y=bquote(f[RB-Harris-LLoG](x))) -> gg1
print(gg1)

ggsave("Plots/LLoG_PDF.eps", gg1, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution







######### hrts B-Harris-Log-Logistic ###########
hrt=function(x,delta,v,theta,c){
  y=((1/gamma(delta))*(theta^(1/v))*(c*(x^(c-1))*(1+x^c)^(-2))*((-log(1-((theta*(1+x^c)^(-v))/(1-(1-theta)*(1+x^c)^(-v)))^(1/v)))^(delta-1))/((1-(1-theta)*(1+x^c)^(-v)))^(1+(1/v))
     
  )/(pgamma(-log(1-((theta*(1+x^c)^(-v))/(1-(1-theta)*(1+x^c)^(-v)))^(1/v)),delta)
  )
  return(y)
}

x_sim=seq(0.01,1.5,by=0.001)
hratios = data.frame(x=x_sim,
                       y1=hrt(x_sim,0.3,0.1,0.9,0.1),
                       y2=hrt(x_sim,2.5,0.7,1.4,1.5),
                       y3=hrt(x_sim,4.5,4.0,19.0,1.3),
                       y4=hrt(x_sim,2.8,9.9,9.8,1.0),
                       y5=hrt(x_sim,1.4,0.6,0.8,5.8)) %>%
  pivot_longer(cols = -x,names_to = "Hratios", values_to = "y")



hratios %>% 
  ggplot(aes(x=x,y=y,colour=Hratios)) +
  geom_line(lwd=.8) +
  scale_colour_manual("",
                      values = c("red", "magenta", "green","cyan","black","blue"),
                      labels = c(expression(delta ~ "=0.3, v =0.1,"~ theta ~"=0.9, c =0.1"),
                                 expression(delta ~ "=2.5, v =0.7,"~ theta ~"=1.4, c =1.5"),
                                 expression(delta ~ "=4.5, v =4.0,"~ theta ~"=19.0, c =1.3"),
                                 expression(delta ~ "=2.8, v =9.9,"~ theta ~"=9.8, c =1.0"),
                                 expression(delta ~ "=1.4, v =0.6,"~ theta ~"=0.8, c =5.8")))+
  scale_linetype_manual("",
                        values = c(1:5),
                        labels = c(expression(delta ~ "=0.3, v =0.1,"~ theta ~"=0.9, c =0.1"),
                                   expression(delta ~ "=2.5, v =0.7,"~ theta ~"=1.4, c =1.5"),
                                   expression(delta ~ "=4.5, v =4.0,"~ theta ~"=19.0, c =1.3"),
                                   expression(delta ~ "=2.8, v =9.9,"~ theta ~"=9.8, c =1.0"),
                                   expression(delta ~ "=1.4, v =0.6,"~ theta ~"=0.8, c =5.8"))) +
  ylim(0,4.5)+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 15,color="black"), #axis font 
        axis.title = element_text(size=15,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        #legend.title = element_text(face="bold",size = 15),
        legend.text.align = 0,
        legend.position = c(0.025,1.0), #position of the legend
        legend.justification = c("left", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        #legend.spacing.x = unit(1.0, 'cm'),
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) + #make transparent background in legend
  labs(y=bquote(h[RB-Harris-LLoG](x))) -> gg2
print(gg2)

ggsave("Plots/LLoG_HRF.eps", gg2, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution


#library(ggpubr)
#gg= ggarrange(gg1,gg2,
#              ncol = 1,
#          labels = "AUTO"
#            )
#print(gg)
#ggsave("pdfhrf_LLoG.eps", gg, width=18, height=25, units="cm", dpi=1080) #save with a specific dimension and resolution

