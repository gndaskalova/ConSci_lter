# Is there a change in beta diversity across the four sample sites?
# Is there a change in beta diversity over time?
# 07/11/17, Conservation Science Exercise, Clara, Jack, Laura, Gergana

# Libraries ----
library(dplyr)
library(ggplot2)
library(vegan)
library(gridExtra)

# Define ggplot2 function
theme_consci <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          axis.title.x=element_text(size=16, face="plain"),             
          axis.title.y=element_text(size=16, face="plain"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=20, vjust=1, hjust=0.5),
          legend.text = element_text(size=12, face="italic"),          
          legend.title = element_text(size=12, face="bold"),
          legend.background = element_rect(color = "black", fill = "transparent", 
                                           size = 4, linetype="blank"))
}

# Load data ----
lter <- read.csv("lter.csv")

# Data formatting ----
# Function to turn data into a presence/absence matrix
splist2presabs <- function(data, sites.col, sp.col, keep.n = FALSE) {
  # version 1.1 (7 May 2013)
  # data: a matrix or data frame with your localities and species (each in a different column)
  # sites.col: the name or index number of the column containing the localities
  # sp.col: the name or index number of the column containing the species names or codes
  # keep.n: logical, whether to get in the resulting table the number of times each species appears in each locality; if false (the default), only the presence (1) or absence (0) are recorded
  
  stopifnot(
    length(sites.col) == 1,
    length(sp.col) == 1,
    sites.col != sp.col,
    sites.col %in% 1 : ncol(data) | sites.col %in% names(data),
    sp.col %in% 1 : ncol(data) | sp.col %in% names(data),
    is.logical(keep.n)
  )
  
  presabs <- table(data[ , c(sites.col, sp.col)])
  presabs <- as.data.frame(unclass(presabs))
  if (!keep.n)  presabs[presabs > 1] <- 1
  presabs <- data.frame(row.names(presabs), presabs)
  names(presabs)[1] <- names(subset(data, select = sites.col))
  rownames(presabs) <- NULL
  return(presabs)
} 

# Pair-wise similarity analysis ----
# Similarity analysis for 1981
lter1981 <- filter(lter, Year == "1981")
lter1981 <- splist2presabs(lter1981, sites.col = 1, sp.col = 3)

sites <- c("Site A1", "Site B1", "Site C1", "Site D1")
rownames(lter1981) <- sites
lter1981 <- lter1981[,-1]

# Calculating Jaccard index of similarity, 0 completely different, 1 the same species
vegdist(lter1981, method="jaccard")

# Similarity analysis for 1996
lter1996 <- filter(lter, Year == "1996")
lter1996 <- splist2presabs(lter1996, sites.col = 1, sp.col = 3)

rownames(lter1996) <- sites
lter1996 <- lter1996[,-1]

# Calculating Jaccard index of similarity, 0 completely different, 1 the same species
vegdist(lter1996, method="jaccard")

# Similarity analysis for 2013
lter2013 <- filter(lter, Year == "2013")
lter2013 <- splist2presabs(lter2013, sites.col = 1, sp.col = 3)

rownames(lter2013) <- sites
lter2013 <- lter2013[,-1]

# Calculating Jaccard index of similarity, 0 completely different, 1 the same species
vegdist(lter2013, method="jaccard")


# Turnover for each site through time ----

# Similarity analysis for A1
A1 <- filter(lter, Site == "Site A1")
A1 <- splist2presabs(A1, sites.col = 2, sp.col = 3)

years <- c("1981", "1996", "2013")
rownames(A1) <- years
A1 <- A1[,-1]

# Calculating Jaccard index of similarity, 0 completely different, 1 the same species
vegdist(A1, method="jaccard")


# Similarity analysis for B1
B1 <- filter(lter, Site == "Site B1")
B1 <- splist2presabs(B1, sites.col = 2, sp.col = 3)

years <- c("1981", "1996", "2013")
rownames(B1) <- years
B1 <- B1[,-1]

# Calculating Jaccard index of similarity, 0 completely different, 1 the same species
vegdist(B1, method="jaccard")

# Similarity analysis for C1
C1 <- filter(lter, Site == "Site C1")
C1 <- splist2presabs(C1, sites.col = 2, sp.col = 3)

years <- c("1981", "1996", "2013")
rownames(C1) <- years
C1 <- C1[,-1]

# Calculating Jaccard index of similarity, 0 completely different, 1 the same species
vegdist(C1, method="jaccard")

# Similarity analysis for D1
D1 <- filter(lter, Site == "Site D1")
D1 <- splist2presabs(D1, sites.col = 2, sp.col = 3)

years <- c("1981", "1996", "2013")
rownames(D1) <- years
D1 <- D1[,-1]

# Calculating Jaccard index of similarity, 0 completely different, 1 the same species
vegdist(D1, method="jaccard")

# Visualising similarity across years for pairs of sites ----
jaccard <- read.csv("jaccard.csv")
jaccard.site <- read.csv("jaccard2.csv")

(similarity <- ggplot(jaccard, aes(x = year, y = jaccard, colour = `Site pair`)) +
    geom_point(size = 4.5, alpha = 0.8) + 
    geom_line(size = 1) +
    geom_point(shape = 1, size = 4.5, colour = "black") +
    labs(x = "Year", y = "Jaccard index of similarity", title = "a) Community similarity") +
    scale_x_continuous(breaks = c(1981, 1996, 2013)) +
    theme_consci())

(similarity2 <- ggplot(jaccard.site, aes(x = Site, y = Jaccard, colour = Period)) +
    geom_point(size = 4.5, alpha = 0.8) + 
    geom_line(size = 1) +
    geom_point(shape = 1, size = 4.5, colour = "black") +
    labs(x = "Year", y = "Jaccard index of similarity", title = "b) Turnover") +
    #scale_x_continuous(breaks = c(1981, 1996, 2013)) +
    theme_consci())

panel <- grid.arrange(similarity, similarity2, ncol = 2)
ggsave(panel, file = "lter_sim_panel.png", width = 10, height = 5)
