
## Breast data
library(devtools)
load_all()

alignNetwork(expression_file = "./data/Breast_expression.csv", adjacency_file = "./data/Breast_network.csv",tissue = "breast")


alignNetwork(expression_file = "./data/Blood_expression.csv", adjacency_file = "./data/Blood_network.csv",tissue = "blood")
