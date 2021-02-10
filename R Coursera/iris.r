
library("dplyr")
library("ggplot2")
library("magrittr")

# display iris dataset
p <- iris
print(p)

# calculate the average sepal width by species
iris %>%
 group_by(Species) %>%
 summarise(Sepal.Width.Avg = mean(Sepal.Width)) %>%
 arrange(Sepal.Width.Avg)
 
 # import the visualization R language library ggplot2
 
 
 # plot the Sepal.Width vs. Sepal.Length
ggplot(data=iris, aes(x=Sepal.Length, y=Sepal.Width, color=Species)) + geom_point(size=3)


options(jupyter.plot_mimetypes = "image/png") 
p <- ggplot(data=iris, aes(x=Sepal.Length, y=Sepal.Width, color=Species)) + geom_point(size=3)
p + ggtitle(sprintf("Mime type = '%s'", getOption("jupyter.plot_mimetypes")))

options(jupyter.plot_mimetypes = "image/png") 
plot(iris)


