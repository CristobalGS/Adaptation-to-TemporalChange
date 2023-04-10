library(xlsx)
library(ggplot2)

# handy preview file for ggsave() (from here: https://gist.github.com/tjmahr/1dd36d78ecb3cff10baf01817a56e895)
ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

file <- "~/Lit review table.xlsx"
table <- read.xlsx(file, sheetIndex = 1, rowIndex = 3:58, colIndex = c(1:11, 13:18, 20:24))

# 1) Plot year distribution
hist(table$Year)
gg.yrs <- ggplot(table, aes(x = Year)) +
  geom_histogram(breaks = seq(1960, 2030, by = 10), col = "black", fill = "gray90") +
  labs(y = "Number of studies") +
  scale_x_continuous(breaks = seq(1960, 2030, by = 10)) +
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
gg.yrs

# 2) Trait: single vs. multivariate and plasticity vs. no plasticity
trait.sum <- table(table$Single.Multivariate, table$Plasticity)[2:4, ]

# single <- sum(table$Single.Multivariate == "S")
# single*100/nrow(table) #78% single
# multi <- sum(table$Single.Multivariate == "M")
# multi*100/nrow(table) #11% multi
# fitness <- sum(table$Single.Multivariate == "F")
# fitness*100/nrow(table) #7% fitness
trait.df <- data.frame(Trait = rep(c("Single", "Multivariate", "Fitness"), 2),
                       Plasticity = c(rep("Yes", 3), rep("No", 3)),
                       Count = c(trait.sum[3, 2], trait.sum[2, 2], trait.sum[1, 2], trait.sum[3, 1], trait.sum[2, 1], trait.sum[1, 1]))
# sum(trait.df[1:3, 3])*100/nrow(table) #51% modelled plasticity
trait.df$Trait <- factor(trait.df$Trait, levels = c("Single", "Multivariate", "Fitness"))
trait.df$Plasticity <- factor(trait.df$Plasticity, levels = c("Yes", "No"))
gg.trait <- ggplot(trait.df, aes(x = Trait, y = Count, fill = Plasticity)) +
  geom_bar(position = "stack", stat = "identity", col = "black") +
  scale_fill_manual(values = c("gray50", "white")) +
  labs(x = "Trait", y = "Number of studies") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.8),
        axis.title.y = element_blank())
gg.trait

# 3) Percentage of studies that modelled plasticity as developmental vs. labile
dev <- table(table$Labile.Developmental.Both, table$Plasticity)[2, 2]
lab <- table(table$Labile.Developmental.Both, table$Plasticity)[3, 2]
both <- table(table$Labile.Developmental.Both, table$Plasticity)[1, 2]
dev*100/(sum(table$Plasticity)) #75% modelled developmental plasticity
lab*100/(sum(table$Plasticity)) #18% modelled labile plasticity
both*100/(sum(table$Plasticity)) #7% modelled both

# 4) Types of temporal environmental change
counts <- c(sum(table$Abrupt.shift), 
            sum(table$Trend), 
            sum(table$Deterministic.cycle), 
            sum(table$Stochastic.change..with.without.autocorrelation.), 
            sum(table$Other != 0))
env <- data.frame(Type = c("Abrupt shift", "Trend", "Deterministic cycle", "Stochastic", "Other"),
                  Count = counts)
env$Type <- factor(env$Type, levels = c("Abrupt shift", "Trend", "Deterministic cycle", "Stochastic", "Other"))
gg.env <- ggplot(env, aes(x = Type, y = Count)) +
  geom_bar(stat = "identity", col = "black", fill = "gray90") +
  labs(x = "Type of temporal\n environmental change", y = "Number of studies") +
  theme_classic() +
  scale_x_discrete(labels = c("Abrupt shift", "Trend", "Cyclic", "Stochastic", "Other")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_blank())
gg.env

# env[1, 2]*100/nrow(table) #11% abrupt shift
# env[2, 2]*100/nrow(table) #35% trend change
# env[3, 2]*100/nrow(table) #18% cycle
# env[4, 2]*100/nrow(table) #51% stochastic change
# env[5, 2]*100/nrow(table) #31% other

# 5) Studies including more than 1 type of change, plus combined
table(rowSums(table[ , 12:15])) # Distribution of studies with 1, 2, 3, or 4 types of change
sum(table$Combined.changes[which(rowSums(table[ , 12:15]) == 2)]) # Number of studies that include 2 types of change and combines them
sum(table$Combined.changes[which(rowSums(table[ , 12:15]) == 3)]) # Number of studies that include 3 types of change and combines them
changes <- data.frame(No = c("1", "2", "3", "4"),
                      Combined = c(rep("Yes", 4), rep("No", 4)),
                      Count = c(0, 7, 0, 0, 31, 2, 1, 0))
changes$Combined <- factor(changes$Combined, levels = c("Yes", "No"))
gg.changes <- ggplot(changes, aes(x = No, y = Count, fill = Combined)) +
  geom_bar(stat = "identity", position = "stack", col = "black") +
  labs(x = "Number of types of\n environmental change per study", y = "Number of studies") +
  scale_fill_manual(values = c("gray50", "white")) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.8),
        axis.title.y = element_blank())
gg.changes

# sum(changes[changes$No == 1, 3])*100/nrow(table) #56% only one type of change
# sum(changes[changes$No == 2, 3])*100/nrow(table) #16% only two types of change
# sum(changes[changes$No == 3, 3])*100/nrow(table) #2% only two types of change

# changes[changes$No == 2 & changes$Combined == "Yes", 3]*100/sum(changes[changes$No == 2, 3]) #78% of the ones that used 2 types combined them

# 6) Percentage of studies that explore:
# A) Generation time: 2%
sum(table$Generation.time)*100/nrow(table)

# B) Costs of plasticity: 13%
sum(table$Costs.of.plasticity)*100/nrow(table)

# c) Genetic correlations: 13%
sum(table$Genetic.correlations)*100/nrow(table)


## Combine all plots into single figure
library(ggpubr)
figure <- ggarrange(gg.yrs, gg.trait, gg.env, gg.changes, 
                    labels = c("A", "B", "C", "D"), ncol = 4, nrow = 1, align = "h",
                    font.label = list(face = "plain"), vjust = 1.25)
figure

# ggpreview(scale = 1, width = 30, height = 9, units = "cm", limitsize = F)
# ggsave("~/lit_rev_fig.pdf",
#        plot = figure, device = cairo_pdf, scale = 1, width = 30, height = 9, units = "cm", limitsize = F)
# ggsave("~/lit_rev_fig.png",
#        plot = figure, device = png, scale = 1, width = 30, height = 9, units = "cm", limitsize = F)
