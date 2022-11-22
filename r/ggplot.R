#-----------------
## Tutoriel ggplot
#-----------------
# Basé sur les codes de Marie-Pierre Etienne

# Slides: https://marieetienne.github.io/premierspas/_presentation/visu.html#1
# Github repository: https://github.com/MarieEtienne/premierspas/

## Load packages
library(GGally)
library(gganimate)
library(ggplot2)
library(ggpubr)
library(plotly)
library(tidyverse)

# remotes::install_github("allisonhorst/palmerpenguins")
# remotes::install_github("wesanderson")

## Load data
data(penguins,package = "palmerpenguins")

## Renommer les colonnes
penguins <- penguins %>%
  rename(bill_l = bill_length_mm, bill_d = bill_depth_mm, flip_l = flipper_length_mm, bm = body_mass_g)
penguins %>%
  print(n=2)

## Grouper et moyenner
penguins %>%
  group_by(species, sex, year, island) %>%
  mutate(n = n()) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  print(n=10)


## Différents types de graphique
#-------------------------------
## Histogramme
penguins %>%
  ggplot() + 
  aes(x = bill_l) + 
  geom_histogram()


## En fréquence
penguins %>%
  ggplot() + 
  aes(x = bill_l) + 
  geom_histogram() +
  aes(y = ..density..)

# Ou autrement
penguins %>%
  ggplot() + 
  geom_histogram(aes(x= bill_l, y = ..density..)) # <--

## Nuage de point
penguins %>% 
  ggplot() + 
  aes( x= bill_l, y=bill_d) + 
  geom_point()


## Boite à moustache
penguins %>% 
  ggplot() + 
  aes( x= species, y=bill_d) + 
  geom_boxplot()


## Enrichir un dessin progressivement
#------------------------------------
# couleur, annotations, theme, etc.

## Pour un nuage de point
gg <- penguins %>%
  ggplot() + 
  aes(x = bill_l) + 
  geom_histogram() 
gg + 
  aes(y = ..density..)

gg <- penguins %>% 
  ggplot() +
  aes( x= bill_l, y=bill_d) +
  geom_point() 
gg <- gg + 
  aes(col = species)  # Ajout de la couleur

## Pour un histogramme
penguins %>%
  ggplot() + aes(x = bill_l, fill = species) + 
  geom_histogram() 

# Faire 3 histogrammes différents
penguins %>%
  ggplot() + aes(x = bill_l, fill = species) + 
  geom_histogram(position = 'identity') 

# Une autre palette de couleur
gg + 
  scale_color_viridis_d()

# Et d'autres palettes sont disponibles
color_darj <- wesanderson::wes_palette(name = "Darjeeling1")
gg + 
  scale_color_manual(values = color_darj)

penguins %>% ggplot() +
  aes(x= bm, y = flip_l) + 
  geom_point() + 
  aes(col = as_factor(year)) +
  scale_color_manual(values = wesanderson::wes_palette(name = "Darjeeling1")) 

gg <- gg + scale_color_manual(values = color_darj) 

## Annotations des axes
gg  +
  labs( x = 'Bill length in mm') +
  labs(y = 'Bill depth in mm') +
  labs(color = "Species")

## Theme
gg <- gg + scale_color_manual(values = color_darj) +
  labs( x = 'Bill length in mm',  y = 'Bill depth in mm',  color = "Species") +
  theme_light()

# Position/forme de la légende
gg +
  theme(legend.position="bottom")

gg + 
  theme(legend.position=c(.9, .6))

gg + theme(legend.position=c(.9, .6),
           text = element_text(size = 10, face = "italic"),
           axis.text.x = element_text(angle=90, hjust=1),
           legend.text = element_text(size = 9, face = 'plain'),
           legend.title = element_text(size = 11, face = "bold") )

# Fixer un thème par défaut
theme_set(theme_light())
theme_update(legend.position="bottom") #BREAK
penguins %>% 
  ggplot() +
  aes( x= bill_l, y=bill_d, col = species) +
  geom_point() #BREAK


## Et d'autres analyses stats
# Regression
gg + 
  geom_smooth(method = 'lm', se = FALSE) +
  geom_smooth(method = 'loess', se = TRUE)

penguins %>% ggplot() + 
  aes(x= bm, y = flip_l, col = as_factor(year) ) + 
  geom_point() + geom_smooth(method = 'lm') +
  scale_color_manual(values = wesanderson::wes_palette(name = "Darjeeling1"))

# Comparaison de plusieurs variables
penguins %>%
  ggplot() + aes(x = bill_l, fill = species) + 
  geom_histogram() +
  facet_wrap(~species) + 
  labs( x = 'Bill length in mm') +
  scale_fill_manual(values = wesanderson::wes_palette('Darjeeling1')) 

# Ajout courbe de densité
penguins %>%  ggplot() + aes(x = bill_l, y = ..density..) +
  facet_wrap(~species) + 
  geom_histogram(alpha=0.5, aes( fill = species)) +
  geom_density(aes(col = species)) + 
  labs( x = 'Bill length in mm') + #BREAK
  scale_fill_manual(values = wesanderson::wes_palette('Darjeeling1')) +  #BREAK
  scale_color_manual(values = wesanderson::wes_palette('Darjeeling1'))   

# Violin plot
penguins %>%  ggplot() + aes(x = species,  y = bill_l) +
  geom_violin(alpha=0.5, aes( fill = species)) +
  labs( y = 'Bill length in mm') + #BREAK
  scale_fill_manual(values = wesanderson::wes_palette('Darjeeling1')) 

# Boxplot + points de données
penguins %>%  ggplot() + aes(x = species,  y = bill_l) +
  geom_boxplot(alpha=0.5, aes( fill = species)) +
  geom_jitter(color="black", size=0.4, alpha=0.8) +
  labs( y = 'Bill length in mm') + #BREAK
  scale_fill_manual(values = wesanderson::wes_palette('Darjeeling1')) 

# Symboles dans le titres des axes (exposants, etc.)
penguins %>%  mutate(mu = bill_l * bill_d) %>% 
  ggplot() + aes(y= mu ) +
  geom_boxplot(alpha=0.5, aes( fill = species)) +
  labs( y = bquote(mu~(mm^2))) +
  scale_fill_manual(values = wesanderson::wes_palette('Darjeeling1')) 


## Package intéressant pour publication

## ggubr
gg_p1 <- gg

gg_p2 <- penguins %>%
  ggplot()  + aes(x = bill_l, y = ..density..) + geom_histogram(alpha=0.5, aes( fill = species)) +
  geom_density(aes(col = species)) +
  labs( x = 'Bill length in mm') +
  scale_fill_manual(values = wesanderson::wes_palette('Darjeeling1')) +
  scale_color_manual(values = wesanderson::wes_palette('Darjeeling1'))

ggpubr::ggarrange(gg_p1, gg_p2, nrow=2, ncol = 1)

ggpubr::ggarrange(gg_p1, gg_p2, nrow=2, ncol = 1, common.legend = TRUE)

## GGally
penguins %>% 
  ggpairs(columns = c(1,3,4,5), mapping = aes(col = species)) +
  scale_color_manual(values = wesanderson::wes_palette('Darjeeling1'))+
  scale_fill_manual(values = wesanderson::wes_palette('Darjeeling1')) + theme(text = element_text(size = 6))

# Control size
penguins %>%
  ggpairs(columns = c(1,3,4,5), mapping = aes(col = species),
                     upper = list(continuous = wrap( "cor",size = 2)),
                     lower = list(continuous = wrap('points', size = .5))) +
  scale_color_manual(values = wesanderson::wes_palette('Darjeeling1'))+
  scale_fill_manual(values = wesanderson::wes_palette('Darjeeling1')) + theme(text = element_text(size = 6))

## gganimate
gg  +
  transition_states(year)  +
  geom_text(x = 56 , y = 15,
            aes(label = as.character(year)),
            size = 8, col = "grey50") +
  theme(legend.position="bottom") 

## plotly
gg  %>%  ggplotly()
