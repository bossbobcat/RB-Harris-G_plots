f1=function(x,delta,v,theta){
  y=(1/gamma(delta))*(-log(1-((theta*(1-pnorm(x,mean(x),sd(x)))^v)/((1-(1-theta)*(1-pnorm(x,mean(x),sd(x)))^v))^(1/v))))^(delta-1)*((theta^(1/v)*dnorm(x,mean(x),sd(x)))/(1-((1-theta)*(1-pnorm(x,mean(x),sd(x)))^v))^(1+(1/v)))
  return(y)
}
x=seq(0,1,by=0.001)

#y1=f1(x, 3.0,.4,1.9)
#y2=f1(x,.8,4.0,1.9)
#y3=f1(x,0.1,3.0,.4)
#y4=f1(x,.8,5.0,1.2)
#y5=f1(x,2.0,1.6,0.4)
#M=data.frame(x,y1,y2,y3,y4,y5)

y1=dnorm(x,mean(x),sd(x))
M=data.frame(x,y1)

library(patchwork)
library(ggplot2)
library(tidyverse)
M %>%
  drop_na() %>%
  pivot_longer(cols = -x,names_to = "para",values_to = "density") ->M #save x and y into M in a systemmatic way

gg1<-ggplot(M, aes(x=x,y=density,col=para)) +     #plot all y's at once
  geom_line(lwd=1.4) +
  scale_color_discrete(name=NULL, #legend name with a title "Parameters" if you want
                       labels = c(expression(paste(mean,'=0.5, ',sd,'=0.2891'))))+ 
  ylim(0,2)+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 30,color="black"), #axis font 
        axis.title = element_text(size=30,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 18),  #legend font
        legend.title = element_text(face="bold",size = 19))+
  theme(
    legend.position = c(1, 1.02), #position of the legend
    legend.justification = c("right", "top"), # location of the legend within the specified location
    legend.box.just = "right",
    legend.margin = margin(11.25, 1.25, 1.2, 1.3) #change this will change the layout of the legend in the box created
  )+
  theme(axis.line.x.bottom = element_line(color = "black", #control axis and its color  
                                          size = 0.8),
        axis.line.y.left = element_line(color = "black", #control axis
                                        size = 0.8))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ #remove all grids, delete if you want grid mesh on
  theme(legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg1)


#ggsave("unidens.eps", gg1, width=18, height=17, units="cm", dpi=1080) #save with a specific dimension and resolution

