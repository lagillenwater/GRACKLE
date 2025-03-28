# GRACKLE

Here, we present a novel NMF approach that applies Graph Regulari-zation Across Contextual KnowLedgE (GRACKLE) to learn latent representations. GRACKLE constrains the unsupervised NMF model using both sample similarity and gene similarity matrices. The model is flexible to incorporate multiple forms of sample metadata, such as phenotype or subtype labels, as well as molecular relationships, such as gene-gene interactions or pathway annotation. 

A schematic of GRACKLE is shown below. 

<img width="400" alt="image" src="https://github.com/user-attachments/assets/bd3436bb-644b-486b-8f52-034d81ab54ff" />

Before matrix decomposition, the sample metadata is transformed into the sample similarity matrix $$S_S$$ by taking the dot product. Similarly, the gene interactions constitute the gene similarity matrix $$S_G$$. These similarity matrices are used to calculate the graph Laplacians $$L_S$$ and $$L_G$$, which are incorporated into the loss function as,

$$\min_{W,H} ||Y - WH||_F^2 + \lambda_1 \text{TR}(W^T L_S W) + \lambda_2 \text{TR}(H^T L_G H) \quad $$

where $$\lambda_1$$ and $$\lambda_2$$ are tuning parameters for graph regularization. Accordingly, the multiplicative update functions from Equation (4) are changed to,

$$w_{ik} \leftarrow w_{ik} \times \frac{YH^T + \lambda_1 S_S W}{WHH^T + \lambda_1 D_S W}, \quad h_{jk} \leftarrow h_{jk} \times \frac{W^T Y + \lambda_2 S_G H}{HWW^T + \lambda_2 D_G H} \quad $$

The model is iteratively trained until a stopping criterion of relative change in $$H$$ is less than $$1 \times 10^{-4}$$ between iterations or a maximum of 100 iterations.
