#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         UNIT PIPE PNEUMATIC MODEL 
#             Melvin T. Tyree       
#       R version: Luciano Pereira & Gabriel Pires
#               09/11/2020          
#     corrected 09/06/2021
#     modified  20/07/2021
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(tidyr)
library(ggplot2)
library(magick)

source("upp.R")

R = 8.3144621/1000    #gas constant
temp = 298.15         #temperature kelvins (25?C)

#Esse loop fiz com voc? no meet
data_all <- data.frame()
for (j in seq(10)) {
  for (i in c(0.05, 0.005,0.0005)) {
    data <- upp(dt = i, PLC = j/10, total_time = 15) %>%
      select(dt1, Cg1) %>%
      mutate(dt_loop = i, plc_loop = j)
    data_all <- rbind(data_all, data)
  }
}

ggplot(data_all) +
  geom_point(aes(dt1, Cg1, col = plc_loop))+
  facet_wrap(~dt_loop)

data <- upp(rings = 6, intact = 12, dt = 0.5)

data_cg <- data %>%
  dplyr::select(dt1, Cg1, Cg2) %>% 
  tidyr::pivot_longer(
    cols = c(Cg1, Cg2),
    names_to = "vessel",
    names_prefix = "gc",
    values_to = "concentration"
  ) %>% 
  mutate(PLC = 0.1)



data_rings <- list()
data_rings_plot <- data.frame()
time_rings <- vector()

for (i in seq(20)) {
  
  data <- upp(rings = i + 2)

  data_cg <- data %>%
    dplyr::select(dt1, starts_with("Cg")) %>% 
    tidyr::pivot_longer(
      cols = starts_with("Cg"),
      names_to = "vessel",
      names_prefix = "gc",
      values_to = "value"
    ) %>% 
    mutate(rings = i)
  
  data_rings_plot <- rbind(data_rings_plot, data_cg)
  data_rings[[i]] <- data
}

qplot(seq_along(time_rings), time_rings)

for (i in 1:20) {
  p<-ggplot(data = data_rings_plot %>%  filter(rings == i),
            aes(x = dt1, y = value*R*temp, col = vessel))+
    geom_point(size=1) +
    scale_x_continuous("Time (s)") +
    scale_y_continuous ("Pressure (kPa)")+
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),
                       aspect.ratio = 1) + # formato final
    theme(text = element_text(size=20),
          axis.title.x = element_text(vjust = - 1.0),
          axis.title.y = element_text(vjust = + 1.0)) +
    labs(title = paste0("Number of rings: ", i + 2))
  ggsave(paste0("rings/",str_pad(i, 3, pad = "0"),"_rings.png"))
}

## list file names and read in
imgs <- list.files("./rings", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 5)

## view animated image
#img_animated

## save to disk
image_write(image = img_animated,
            path = "rings.gif")

data_intact <- list()
data_intact_plot <- data.frame()
time_intact <- vector()

for (i in seq(20)) {
  print(paste0("Starting ",i," of 20"))
  start.time <- Sys.time()
  data <- upp(intact = i + 2)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_intact[i] <- time.taken
  print(paste0("Done in ", round(time.taken,2), " sec"))
  
  data_cg <- data %>%
    dplyr::select(dt1, starts_with("v1")) %>% 
    tidyr::pivot_longer(
      cols = starts_with("v1"),
      names_to = "vessel",
      names_prefix = "v1",
      values_to = "value"
    ) %>% 
    mutate(intact = i)
  
  data_intact_plot <- rbind(data_intact_plot, data_cg)
  data_intact[[i]] <- data
}

qplot(seq_along(time_intact), time_intact)

for (i in 1:20) {
  p<-ggplot(data = data_intact_plot %>%  filter(intact == i),
            aes(x = dt1, y = value, col = vessel))+
    geom_point(size=1) +
    scale_x_continuous("Time (s)") +
    scale_y_continuous (bquote('C'[s]*'(mol'~ m^-3*')'))+
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),
                       aspect.ratio = 1) + # formato final
    theme(text = element_text(size=20),
          axis.title.x = element_text(vjust = - 1.0),
          axis.title.y = element_text(vjust = + 1.0)) +
    labs(title = paste0("Number of intact vessel: ", i + 2))
  ggsave(paste0("intact/",str_pad(i, 3, pad = "0"),"_intact.png"))
}

## list file names and read in
imgs <- list.files("./intact", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 5)

## view animated image
#img_animated

## save to disk
image_write(image = img_animated,
            path = "intact.gif")












data_PLC <- list()
data_PLC_plot <- data.frame()
time_PLC <- vector()

for (i in seq(100)) {
  print(paste0("Starting ",i))
  start.time <- Sys.time()
  data <- upp(PLC = i/100)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_PLC[i] <- time.taken
  print(paste0("Done in ", round(time.taken,2), " sec"))
  
  data_cg <- data %>%
    dplyr::select(dt1, starts_with("v1")) %>% 
    tidyr::pivot_longer(
      cols = starts_with("v1"),
      names_to = "vessel",
      names_prefix = "v1",
      values_to = "value"
    ) %>% 
    mutate(PLC = i)
  
  data_PLC_plot <- rbind(data_PLC_plot, data_cg)
  data_PLC[[i]] <- data
}

qplot(seq_along(time_PLC), time_PLC)

for (i in 1:100) {
  p<-ggplot(data = data_PLC_plot %>%  filter(PLC == i),
            aes(x = dt1, y = value, col = vessel))+
    geom_point(size=1) +
    scale_x_continuous("Time (s)") +
    scale_y_continuous (bquote('C'[s]*'(mol'~ m^-3*')'))+
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),
                       aspect.ratio = 1) + # formato final
    theme(text = element_text(size=20),
          axis.title.x = element_text(vjust = - 1.0),
          axis.title.y = element_text(vjust = + 1.0)) +
    labs(title = paste0("PLC: ", i/100))
  ggsave(paste0("PLC/",str_pad(i, 3, pad = "0"),"_PLC.png"))
}

## list file names and read in
imgs <- list.files("./PLC", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 5)

## view animated image
#img_animated

## save to disk
image_write(image = img_animated,
            path = "PLC.gif")







