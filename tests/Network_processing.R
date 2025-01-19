## Breast data
library(devtools)
load_all()
#processNetwork(expression_file = "../GRACKLE_data/data/Breast_expression.csv", adjacency_file = "../GRACKLE_data/data/Breast_network.csv",tissue = "breast")

processNetwork(expression_file = "../GRACKLE_data/data/Blood_expression.csv", adjacency_file = "../GRACKLE_data/data/Blood_network.csv",tissue = "blood")
