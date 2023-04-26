library(ggplot2)
library(Cairo)

cat <- read.delim("cat.csv", header=FALSE)
met <- read.table("met.csv", quote="\"", comment.char="")
perf <- read.csv("perf.csv", sep="")

df = data.frame(rate = perf, metric = met, parameters = cat)
colnames(df) = c("Rate", "Metric", "Parameters")
p = ggplot(df, aes(fill=Metric, y=Rate, x=Parameters)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_brewer() + 
  theme_minimal()

ggsave(filename="../Image/perf.pdf", plot=p, device=cairo_pdf,
       width = 8, height = 4, units = "in", dpi = 1200)
