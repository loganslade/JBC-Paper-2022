#####Used for creating Figs. 9, 10, 11, Supp. Figs. 5-8###### 

setwd(dir)

library(ggplot2)
library(dplyr)
library(tibble)
library(forcats)
library(ggridges)
library(readr)

#####Import and filter data from cell profiler results######
hci.edu <- read_csv(file.choose())
head(hci.edu, n = 5)


hci.edu.e <- filter(hci.edu, !treat %in% c("nedu","nab"))

max.e <- max(hci.edu.e$mean_edu)
min <- min(hci.edu.e$mean_edu)

hci.edu.e$treat <- as.factor(hci.edu.e$treat)
levels(hci.edu.e$treat) <- c("siCTRL", "siTFEB#2")
hci.edu.e$exptreat <- interaction(hci.edu.e$treat, hci.edu.e$exp)
hci.edu.e$scenetreat <- interaction(hci.edu.e$treat, hci.edu.e$exp, hci.edu.e$scene)

#######Normalization of intensity value across experiments######
corany <- function(a,b,c){
  corfac <- a  %>% group_by(.data[[c]]) %>% summarise(mean = mean(.data[[b]])) %>%
    mutate(corfac = mean[1]/mean) %>% select(-mean)  
  
  cor <- left_join(x = a, y = corfac) %>% mutate("cor_{b}" := (.data[[b]]*corfac)) %>% select(-corfac)                   
  return(cor)
}

hci.edu.a <- hci.edu.e %>% corany("int_dna", "scenetreat")


######Downsample######
hci.edu.count <- hci.edu.a %>% group_by(exptreat) %>% summarise(n())
hci.edu.count
hci.edu.ds <- hci.edu.a %>% group_by(exptreat) %>% sample_n(2000)
hci.edu.c <- hci.edu.a


#######Univariate density plots######
qc <- ggplot(hci.edu.c, aes(factor(exptreat),cor_int_dna))
qc + geom_boxplot() +
  scale_x_discrete()


dna.dense <- ggplot(hci.edu.c, aes(cor_int_dna, fill = treat)) 
dna.dense + geom_density(adjust=1.1, alpha=0.3) + 
  theme_classic()+
  labs(title="MCF10", x="DNA intensity", y="Density", fill = "siRNA") +  
  scale_x_continuous(limits = c(0, 800))+
  scale_fill_manual( values = c("turquoise3", "red1"))+
  theme(text = element_text(size=16),
        plot.title = element_text(hjust = 0.5))

edu.dense <- ggplot(hci.edu.a, aes(mean_edu, fill = treat)) 
edu.dense + geom_density(adjust=1.2, alpha=0.3) + 
  theme_classic()+
  labs(title="MCF10A", x="MFI EdU", y="Density", fill = "siRNA") +  
  scale_x_continuous(trans = "log10", limits = c(min, max.e))+
  scale_fill_manual( values = c("turquoise3", "red1"))+
  theme(text = element_text(size=16),
        plot.title = element_text(hjust = 0.5))

mcm.dense <- ggplot(hci.edu.c, aes(mean_mcm2, fill = treat)) 
mcm.dense + geom_density(adjust=1.2, alpha=0.3) + 
  theme_classic()+
  labs(title="MCF10A", x="MFI MCM2", y="Density", fill = "siRNA") +  
  scale_x_continuous(limits = c(0, 10000))+
  scale_fill_manual( values = c("turquoise3", "red1"))+
  theme(text = element_text(size=16),
        plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 3000)

######2d and 3d plots#####
gg <- hci.edu.c %>% filter(area < 2500) %>% ggplot(aes(x=cor_int_dna, y=mean_edu))
gg + geom_jitter(aes(colour = treat), alpha=0.4, size = 0.5)+
  scale_y_continuous(trans = "log10", limits = c(min, max.e))+
  scale_colour_manual( values = c("turquoise3", "red1"))+
  xlim(0, 1000)+
  labs(title = "MDA-MB-231", x="DNA Integrated Intensity", y="MFI EdU", colour = "siRNA")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1))+
  facet_grid(cols = vars(treat))

h2ax <- hci.edu.ds %>% filter(area < 2500) %>% ggplot(aes(x=cor_int_dna, y=mean_edu))
h2ax + geom_jitter(aes(fill = mean_mcm2), color = "black", pch = 21, alpha=1, size = 1.9)+
  scale_y_continuous(trans = "log10", limits = c(min*0.999, max.e))+
  scale_fill_gradient2(low = "yellow", mid = "orange", high = "red", na.value = "red", 
                           midpoint = 0.068, limits = c(0,0.1))+
  xlim(0, 1000)+
  labs(title = "MDA-MB-231", x="DNA Integrated Intensity", y="MFI EdU", fill = expression(paste("MFI MCM2")))+
  theme_classic()+
  facet_grid(cols = vars(treat))+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1))

md <- hci.edu.ds %>% filter(area < 2500) %>% ggplot(aes(x=cor_int_dna, y=mean_mcm2))
md +  geom_jitter(aes(colour = treat), alpha=0.1, size = 0.1)+
  scale_x_continuous(limits = c(0, 1000))+
  scale_y_continuous(trans = "log10", limits = c(0.01, 0.4))+
  scale_colour_manual( values =c("turquoise3", "red1"))+
  labs(title = "MDA-MB-231", x="DNA Integrated Intensity", y= expression(paste("MFI MCM2")), fill = expression(paste("MFI ", gamma, "H2AX")))+
  theme_classic()+
  facet_grid(cols = vars(treat))+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1))

########CC phase intensity values##############
x1 <- 150
x2 <- 500
x3 <- x2+(x2-x1)
y1 <- min*0.9999
y2<- 0.006
y3 <- 0.03
y4 <- max.e

#######Gate validation###########

gg.boxes <- hci.edu.c %>% ggplot(aes(x=cor_int_dna, y=mean_edu))
gg.boxes + geom_jitter(aes(colour = "red4"), alpha=0.5, size = 0.5)+
  scale_y_continuous(trans = "log10", limits = c(min*0.9999, max.e))+
  xlim(0, 1000)+
  labs(title = "hci-MB-231", x="DNA Integrated Intensity", y="MFI EdU")+
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "grey88"))+
   annotate("rect", xmin= x1, xmax= x2, ymin= y1, ymax= y2, 
           alpha=0, color = "black", linetype = "dotted")+
   annotate("rect", xmin= x1, xmax= x2, ymin= y2, ymax= y3, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin= x1, xmax= x2, ymin= y3, ymax= y4, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin= x2, xmax= x3, ymin= y3, ymax= y4, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin=x2, xmax=x3, ymin=y2, ymax=y3, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin=x2, xmax=x3, ymin=y1, ymax=y2, 
           alpha=0, color = "black", linetype = "dotted")+ 
  facet_grid(cols = vars(treat), rows = vars(exp))


#####Gates to cell cycle phase#####
hci.edu.tib <- hci.edu.a %>% as_tibble 
gated.hci.edu <- hci.edu.tib %>% mutate(phase = ifelse(between(hci.edu.tib$mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x1,x2), "G1",
                                                             ifelse(between(hci.edu.tib$mean_edu,y2,y4) & between(hci.edu.tib$cor_int_dna,x1,x2), "S1",
                                                                    ifelse(between(hci.edu.tib$mean_edu,y2,y4) & between(hci.edu.tib$cor_int_dna,x2,x3), "S2",
                                                                           ifelse(between(hci.edu.tib$mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x2,x3), "G2", "OT"))))) 


tgated.hci.edu <- gated.hci.edu %>% mutate(phase2 = ifelse(between(hci.edu.tib$mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x1,x2), "G1",
                                                              ifelse(between(hci.edu.tib$mean_edu,y2,y4) & between(hci.edu.tib$cor_int_dna,x1,x3), "S",
                                                                    ifelse(between(hci.edu.tib$mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x2,x3), "G2", "OT"))))


sixgated.hci.edu <- tgated.hci.edu %>% mutate(phase3 = ifelse(between(hci.edu.tib$mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x1,x2), "G1",
                                                        ifelse(between(hci.edu.tib$mean_edu,y2,y3) & between(hci.edu.tib$cor_int_dna,x1,x2), "ES",
                                                            ifelse(between(hci.edu.tib$mean_edu,y3,y4) & between(hci.edu.tib$cor_int_dna,x1,x2), "S1",
                                                                ifelse(between(hci.edu.tib$mean_edu,y3,y4) & between(hci.edu.tib$cor_int_dna,x2,x3), "S2",
                                                                     ifelse(between(hci.edu.tib$mean_edu,y2,y3) & between(hci.edu.tib$cor_int_dna,x2,x3), "LS",
                                                                        ifelse(between(hci.edu.tib$mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x2,x3), "G2", "OT"))))))) 


sixgated.hci.edu <- sixgated.hci.edu %>%filter (phase != "OT") %>%filter(phase2 != "OT") %>% filter(phase3 != "OT")


write.csv(sixgated.hci.edu, "6pGated_MDAsiTFEDU_MCM2_v2.csv")

gated.edu.f <- sixgated.hci.edu
gated.edu.f$phase3 <- factor(gated.edu.f$phase3, 
                            levels = c("G1", "ES", "S1", "S2", "LS", "G2"))

gated.edu.f$treatphase <- interaction(gated.edu.f$treat, gated.edu.f$phase3)

#######Gate validation and intensity distribution by phase######

dg <- gated.edu.f %>% ggplot(aes(mean_mcm2, fill = treat))
dg  + geom_density(adjust=1.2, alpha=0.3) + 
  theme_classic()+
  labs(x="MFI MCM2", y="Density", fill = "siRNA") +  
  scale_x_continuous(limits = c(0, 0.4))+
  scale_fill_manual( values = c("turquoise3", "red1"))+
  theme(text = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(0.75, "lines"),
        axis.text = element_text(colour = "black", size = 12))+
  facet_grid(cols = vars(phase3))


mcmg <- gated.edu.f %>% ggplot(aes(x=cor_int_dna, y=mean_edu))
mcmg + geom_jitter(aes(colour = phase3), alpha=0.8, size = 0.5)+
  scale_y_continuous(trans = "log10", limits = c(min*0.999, max.e))+
  scale_colour_manual( values = c("green", "gold", "red", "violet", "orange", "blue"))+
  xlim(0, 1000)+
  labs(title = NULL, x="DNA Integrated Intensity", y="MFI EdU", colour = "Gate")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1))+
  annotate("rect", xmin= x1, xmax= x2, ymin= y1, ymax= y2, 
         alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin= x1, xmax= x2, ymin= y2, ymax= y3, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin= x1, xmax= x2, ymin= y3, ymax= y4, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin= x2, xmax= x3, ymin= y3, ymax= y4, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin=x2, xmax=x3, ymin=y2, ymax=y3, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin=x2, xmax=x3, ymin=y1, ymax=y2, 
           alpha=0, color = "black", linetype = "dotted")+
  guides(colour = guide_legend(override.aes = list(size=5)))






