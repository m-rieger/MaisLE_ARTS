#### plot nextmatch with dummy data ####
library(ggplot2)

df.filter <- data.frame(signal = c("s0", "s1", "s2", NA, NA, NA),
                        nm = c(NA, "(s0 nextmatch)", "(s1 nextmatch)", NA, NA, NA),
                        y = 150.100,
                        x = c(0, 1.1, 2.15, 0.7, 1.3, 2.65),
                        color = c("darkblue", "darkorange", "grey30", "grey", "grey", "grey"))
texp <- 1 ## expected time interval
x0 <- df.filter$x[which(df.filter$signal == "s0")]
c0 <- df.filter$color[which(df.filter$signal == "s0")]
x1 <- df.filter$x[which(df.filter$signal == "s1")]
c1 <- df.filter$color[which(df.filter$signal == "s1")]
x2 <- df.filter$x[which(df.filter$signal == "s2")]

g <- ggplot() + 
  ## area for s0
  geom_line(aes(x = x0+texp, y = c(150.050, 150.200)), color = c0, lty = "dashed") +
  geom_rect(aes(xmin = x0+texp-0.5*texp, xmax = x0+texp+0.5*texp, ymin = 150.050, ymax = 150.200),
            fill = c0, alpha = 0.2) +
  geom_text(aes(x = x0+texp, y = 150.205), label = "window to search for s0 nextmatch", color = c0, size = 5) +
  ## area for s1
  geom_line(aes(x = x1+texp, y = c(150.050, 150.200)), color = c1, lty = "dashed") +
  geom_rect(aes(xmin = x1+texp-0.5*texp, xmax = x1+texp+0.5*texp, ymin = 150.050, ymax = 150.200),
            fill = c1, alpha = 0.2) +
  geom_text(aes(x = x1+texp, y = 150.205), label = "window to search for s1 nextmatch", color = c1, size = 5) +
  
  ## lines s0
  geom_line(aes(x = c(x0, x0+texp), y = 150.125), lwd = 1) +
  geom_text(aes(x = mean(c(x0, x0+texp)), y = 150.13), label = "texp", size = 5) +
  
  geom_line(aes(x = c(x0, x1), y = 150.135), color = c0, lwd = 1) +
  geom_text(aes(x = mean(c(x0, x1)), y = 150.14), color = c0, label = "tact", size = 5) +
  
  geom_line(aes(x = c(x1, x0+texp), y = 150.145), color = c0, lwd = 1) +
  geom_text(aes(x = mean(c(x1, x0+texp)), y = 150.15), color = c0, label = "Nextmatch delta", size = 5) +
  
  ## lines s1
  geom_line(aes(x = c(x1, x1+texp), y = 150.125), lwd = 1) +
  geom_text(aes(x = mean(c(x1, x1+texp)), y = 150.13), label = "texp", size = 5) +
  
  geom_line(aes(x = c(x1, x2), y = 150.135), color = c1, lwd = 1) +
  geom_text(aes(x = mean(c(x1, x2)), y = 150.14), color = c1, label = "tact", size = 5) +
  
  geom_line(aes(x = c(x2, x1+texp), y = 150.145), color = c1, lwd = 1) +
  geom_text(aes(x = mean(c(x2, x1+texp)), y = 150.15), color = c1, label = "Nextmatch delta", size = 5) +
  
  geom_point(aes(x = x, y = y), size = 5, color = df.filter$color, data = df.filter) +
  geom_text(aes(label = signal, x = x, y = y-0.005), data = df.filter, size = 5, color = df.filter$color) +
  geom_text(aes(label = nm, x = x, y = y-0.01), data = df.filter[!is.na(df.filter$nm),], size = 5, color = c(c0, c1)) +
  
  xlab("time [s]") +
  ylab("frequency [MHz]") +
  ylim(150.05, 150.21) +
  theme_light(base_size = 15)

ggsave("./plotNextmatch.pdf", plot = g, width = 30, height = 18, device = "pdf", units = "cm")
