#install.packages("ggalluvial")
library("ggalluvial")

setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/dataset/Data_final")
al<-read.csv("alluvial_plot_data.csv", header=T,  stringsAsFactors = FALSE, na.strings = c("", NA))
al<-read.csv("alluvial_plot_data.csv", header=T,  stringsAsFactors = TRUE, na.strings = c("", NA))


al <- read.csv("alluvial-plot/alluvial_plot_data.csv", header=T, stringsAsFactors=T, na.strings=c("", NA))
class(al)
str(al)
al<-as_tibble(al)

al$Level_1<-factor(al$Level_1, al$Level_1, ordered = FALSE)
al$Level_2<-factor(al$Level_2, unique(al$Level_2), ordered = FALSE)

# ids
# al<-add_column(al,
#            Level_1_ids=get_onto_id(al$Level_1_names, ONT), 
#            Level_2_ids=get_onto_id(al$Level_2, ONT), 
#            Level_3_ids=get_onto_id(al$Level_3, ONT) )
# 
# write.csv(al, file="alluvial_plot_data.csv",  row.names = FALSE)
# ##

ggplot(data = al,
       aes(axis1 = Level_1, axis2 = Level_2, axis3 = Level_3, weight = Level_1_nstates))+
  scale_x_discrete(limits = c("Level 1", "Level 2", "Level 3")) +
  geom_alluvium(aes(fill = Level_2), color = "black",   reverse = FALSE,  size=.1, width = 0, 
                 knot.pos = 0, na.rm = TRUE) +
  
  guides(fill = FALSE) +
  geom_stratum( reverse = FALSE, fill = c(F$color, F2$color[1:6] ), color = "black", size=.2 , na.rm = TRUE) + 
  geom_text(stat = "stratum", label.strata = TRUE, reverse = FALSE, angle = 60) +
  #geom_label(stat = "stratum", label.strata = TRUE, angle = 60) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal() +
  geom_flow()+
  theme(plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x=element_blank()
  ) +
  coord_flip()+
  ggtitle("")

##
ggplot(data = al,
       aes(axis1 = Level_1, axis2 = Level_2, axis3 = Level_3, weight = Level_1_nstates))+
  scale_x_discrete(limits = c("Level 1", "Level 2", "Level 3")) +
  geom_alluvium(aes(fill = Level_2), width = 0,  reverse = FALSE) +
  
  guides(fill = FALSE) +
  #geom_flow(aes.flow = "backward") +
  geom_stratum(alpha = 1, reverse = FALSE) + 
  geom_text(stat = "stratum", label.strata = TRUE, reverse = FALSE, angle = 0) +
  #geom_label(stat = "stratum", label.strata = TRUE, angle = 0) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal() +
  #coord_flip()+
  ggtitle("")

##


##
data(vaccinations)
levels(vaccinations$response) <- rev(levels(vaccinations$response))
ggplot(vaccinations,
       aes(x = survey, stratum = response, alluvium = subject,
           weight = freq,
           fill = response, label = response)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none") +
  ggtitle("vaccination survey responses at three points in time")

##


titanic_wide <- data.frame(Titanic)
head(titanic_wide)

ggplot(data = titanic_wide,
       aes(axis1 = Class, axis2 = Sex, axis3 = Age,
           weight = Freq)) +
  scale_x_discrete(limits = c("Class", "Sex", "Age")) +
  geom_alluvium(aes(fill = Survived)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("passengers on the maiden voyage of the Titanic",
          "stratified by demographics and survival")
#############


# one axis
ggplot(as.data.frame(Titanic),
       aes(weight = Freq,
           axis = Class)) +
  geom_lode(aes(fill = Class, alpha = Survived)) +
  scale_x_continuous(breaks = 1, labels = c("Class")) +
  scale_alpha_manual(values = c(.25, .75))

gg <- ggplot(as.data.frame(Titanic),
             aes(weight = Freq,
                 axis1 = Class, axis2 = Sex, axis3 = Age,
                 fill = Survived))
# alluvia and lodes
gg + geom_alluvium() + geom_lode()
# lodes as strata
gg + geom_alluvium() +
  geom_stratum(stat = "alluvium")
#
# full axis width
ggplot(as.data.frame(Titanic),
       aes(weight = Freq,
           axis1 = Class, axis2 = Sex, axis3 = Age, axis4 = Survived)) +
  geom_stratum(width = 1) + geom_text(stat = "stratum", label.strata = TRUE) +
  scale_x_continuous(breaks = 1:4,
                     labels = c("Class", "Sex", "Age", "Survived"))

# use of facets
ggplot(as.data.frame(Titanic),
       aes(weight = Freq,
           axis1 = Class, axis2 = Sex)) +
  geom_flow(aes(fill = Survived)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  lode_rightward(n=3, i=2)+
  scale_x_continuous(breaks = 1:2, labels = c("Class", "Sex")) +
  facet_wrap(~ Age, scales = "free_y")


###################3

dat_ggforce2 <- dat_ggforce %>%
  mutate_at(vars(Modeled, Tested), 
            funs(factor(., levels = fct_levels)))

ggplot(dat_ggforce2, aes(x = x, id = id, split = y, value = freq)) +
  geom_parallel_sets(aes(fill = Tested), alpha = alpha, axis.width = 0.2,
                     n=100, strength = 0.5) +
  geom_parallel_sets_axes(axis.width = 0.25, fill = "gray95",
                          color = "gray80", size = 0.15) +
  geom_parallel_sets_labels(colour = 'gray35', size = 4.5, angle = 0, fontface="bold") +
  scale_fill_manual(values  = c(A_col, C_col, B_col)) +
  scale_color_manual(values = c(A_col, C_col, B_col)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.title.x  = element_blank()
  )


dat_alluvial <- dat %>%
  dplyr::select(Modeled, Tested, freq)

alluvial::alluvial(al[,c(2,3,4)], freq=al$Level_1_nstates,
                   xw=0.2, alpha=1, gap.width=0.2, col=F2$color, cw=.2,
                  
                   blocks=TRUE, cex = 1.1, cex.axis = 1.5
               
)

ggplot(data = al,
       aes(axis1 = Level_1, axis2 = Level_2, axis3 = Level_3, weight = Level_1_nstates))+
  scale_x_discrete(limits = c("Level 1", "Level 2", "Level 3")) +
  geom_alluvium(aes(fill = Level_2),   reverse = FALSE, color = "grey", size=1, width = 0, 
                knot.pos = 0, na.rm = TRUE) +

tit <- as.data.frame(Titanic, stringsAsFactors = FALSE)
head(tit)
alluvial(tit[,1:4], freq=tit$Freq,
         col = ifelse(tit$Survived == "Yes", "orange", "grey"),
         border = ifelse(tit$Survived == "Yes", "orange", "grey"),
         hide = tit$Freq == 0,
         cex = 0.7
)

######333
devtools::install_github('thomasp85/ggforce')
dat_raw <- data.frame(Tested  = sample(c("A","B","C"),100,
                                       replace = TRUE,prob=c(0.2,0.6,0.25)),
                      Modeled = sample(c("A","B","C"),100,
                                       replace = TRUE,prob=c(0.56,0.22,0.85)),
                      stringsAsFactors = FALSE)
dat <- dat_raw %>%
  group_by(Tested,Modeled) %>%
  summarise(freq = n()) %>%
  ungroup()

dat_ggforce <- dat  %>%
  gather_set_data(1:2) %>%        # <- ggforce helper function
  arrange(x,Tested,desc(Modeled))

ggplot(data = al,
       aes(axis1 = Level_1, axis2 = Level_2, axis3 = Level_3, weight = Level_1_nstates))+
  scale_x_discrete(limits = c("Level 1", "Level 2", "Level 3")) +
  geom_alluvium(aes(fill = Level_2), color = "black",   reverse = FALSE,  size=.1, width = 0, 
                knot.pos = 0, na.rm = TRUE) +

ggplot(dat_ggforce, aes(x = x, id = id, split = y, value = freq)) +
  
  geom_parallel_sets(aes(fill = Tested), alpha = alpha, axis.width = 0.2,
                     n=100, strength = 0.5) +
  geom_parallel_sets_axes(axis.width = 0.25, fill = "gray95",
                          color = "gray80", size = 0.15) +
  geom_parallel_sets_labels(colour = 'gray35', size = 4.5, angle = 0, fontface="bold") +
  scale_fill_manual(values  = c(A_col, B_col, C_col)) +
  scale_color_manual(values = c(A_col, B_col, C_col)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.title.x  = element_blank()
  )
