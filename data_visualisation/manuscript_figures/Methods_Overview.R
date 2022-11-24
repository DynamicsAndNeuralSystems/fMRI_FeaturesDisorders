################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(ggseg3d)
library(icesTAF)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/")
icesTAF::mkdir(plot_path)

################################################################################
# Figure 1: methods overview
################################################################################
icesTAF::mkdir(paste0(plot_path, "Figure1/"))

# Brain figure
data.frame(
  region = c("inferior parietal",
             "rostral middle frontal"),
  p = c(0, 1),
  stringsAsFactors = F) %>%
  ggseg3d(colour = "p", text = "p",
          palette = c("deepskyblue2" = 0 ,
                      "darkorange2" = 1),
          show.legend=F) %>% 
  remove_axes() %>%
  pan_camera("left medial")

# TS line plots
set.seed(127)
ar.sim1<-arima.sim(model=list(ar=c(.9,-.2)),n=80)
ar.sim2 <- arima.sim(model=list(ar=c(0.3, 0.1)), n=80)

#Our transformation function
scaleFUN <- function(x) sprintf("%.1f", x)

# Inferior parietal (blue)
data.frame(Time = 1:length(ar.sim1),
           BOLD = as.numeric(ar.sim1),
           Time_Series = "TS1") %>%
  ggplot(data=., mapping=aes(x=Time, y=BOLD, color=Time_Series)) +
  geom_line(size=0.6, alpha=0.9) +
  scale_color_manual(values=c("deepskyblue2")) + 
  scale_y_continuous(labels=scaleFUN) +
  theme(legend.position="none")
ggsave(paste0(plot_path, "Figure1/BOLD_TS_blue.png"), width=2.3, height=1.8, 
       units="in", dpi=1200, bg = "transparent")

# Rostral middle frontal (orange)
data.frame(Time = 1:length(ar.sim2),
           BOLD = as.numeric(ar.sim2),
           Time_Series = "TS2") %>%
  ggplot(data=., mapping=aes(x=Time, y=BOLD, color=Time_Series)) +
  geom_line(size=0.6, alpha=0.9) +
  scale_color_manual(values=c("darkorange2")) +
  scale_y_continuous(labels=scaleFUN) +
  theme(legend.position="none")
ggsave(paste0(plot_path, "Figure1/BOLD_TS_orange.png"), width=2.3, height=1.8, 
       units="in", dpi=1200, bg = "transparent")

# Plots together
data.frame(Time = c(1:length(ar.sim1), 1:length(ar.sim2)),
           BOLD = c(as.numeric(ar.sim1), as.numeric(ar.sim2)),
           Time_Series = c(rep("TS1", length(ar.sim1)),
                           rep("TS2", length(ar.sim2)))) %>%
  ggplot(data=., mapping=aes(x=Time, y=BOLD, color=Time_Series)) +
  geom_line(size=0.6, alpha=0.9) +
  scale_color_manual(values=c("deepskyblue2", "darkorange2")) +
  scale_y_continuous(labels=scaleFUN) +
  theme(legend.position="none")
ggsave(paste0(plot_path, "Figure1/BOLD_TS_both.png"), width=2.3, height=1.8, 
       units="in", dpi=1200, bg = "transparent")

