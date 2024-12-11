library(ggplot2)
library(dplyr)

cond <- expand.grid(t=1:100, 
                    delay = seq(0.2,2, by=0.2),
                    forgetting=seq(0.5,1,by=0.1)) %>% 
  mutate(f = lr_default(t=t, delay=delay, forgetting=forgetting))

head(cond)
ggplot(cond, aes(x=t,y=f, group=forgetting, colour=forgetting))+
  geom_line()+
  facet_grid(delay~.)+
  scale_color_viridis_c()+
  theme_classic()

                