# GRACKLE

Here, we present a novel NMF approach that applies Graph Regulari-zation Across Contextual KnowLedgE (GRACKLE) to learn latent representations. GRACKLE constrains the unsupervised NMF model using both sample similarity and gene similarity matrices. The model is flexible to incorporate multiple forms of sample metadata, such as phenotype or subtype labels, as well as molecular relationships, such as gene-gene interactions or pathway annotation. 

A schematic of GRACKLE is shown below. 

<img width="500" alt="image" src="https://github.com/user-attachments/assets/bd3436bb-644b-486b-8f52-034d81ab54ff" />

Before matrix decomposition, the sample metadata is transformed into the sample similarity matrix $$S_S$$ by taking the dot product. Similarly, the gene interactions constitute the gene similarity matrix $$S_G$$. These similarity matrices are used to calculate the graph Laplacians $$L_S$$ and $$L_G$$, which are incorporated into the loss function as,

$$\min_{W,H} ||Y - WH||_F^2 + \lambda_1 \text{TR}(W^T L_S W) + \lambda_2 \text{TR}(H^T L_G H) \quad $$

where $$\lambda_1$$ and $$\lambda_2$$ are tuning parameters for graph regularization. Accordingly, the multiplicative update functions from Equation (4) are changed to,

$$w_{ik} \leftarrow w_{ik} \times \frac{YH^T + \lambda_1 S_S W}{WHH^T + \lambda_1 D_S W}, \quad h_{jk} \leftarrow h_{jk} \times \frac{W^T Y + \lambda_2 S_G H}{HWW^T + \lambda_2 D_G H} \quad $$

The model is iteratively trained until a stopping criterion of relative change in $$H$$ is less than $$1 \times 10^{-4}$$ between iterations or a maximum of 100 iterations. These model parameters may be fine tuned based on your scientific question.

## Prerequisites

*   **R:**  A statistical computing language and environment.
    *   Download and install the latest version of R from {Link: CRAN https://cran.r-project.org/}

## Download GRACKLE
**Clone the GRACKLE Repository**:
To get started, clone the GRACKLE repository from GitHub:

```bash
git clone https://github.com/your-username/GRACKLE.git
cd GRACKLE
```

## Install Required R Packages

To run GRACKLE, ensure you have the following R packages installed. Follow these steps:

1. **Install Required Packages**:
    Open R and run the following commands to install the necessary packages:
    ```R
    install.packages(c("tidyverse", "igraph", "parallel", "optparse", "devtools", "reticulate", "tensorflow"))
    ```

2. **Set Up Python Environment**:
    GRACKLE uses Python via the `reticulate` package. Ensure you have a virtual environment set up and activate it:
    ```R
    library(reticulate)
    use_virtualenv("/path/to/env")
    ```

3. **Additional Steps To Install TensorFlow**:
    After installing the `tensorflow` R package, run the following command in R to install TensorFlow:
    ```R
    library(tensorflow)
    install_tensorflow()
    ```

4. **Verify Installation**:
    To ensure TensorFlow is installed correctly, run:
    ```R
    library(tensorflow)
    tf$constant("Hello, TensorFlow!")
    ```

5. **GPU Support (Optional)**:
    If you want to enable GPU support, ensure you have the necessary CUDA and cuDNN libraries installed. Then, install TensorFlow with GPU support:
    ```R
    install_tensorflow(version = "gpu")
    ```

For more details, refer to the official TensorFlow for R documentation: [TensorFlow for R](https://tensorflow.rstudio.com/).



## Simulations

To establish a performance benchmark, we compared GRACKLE to other NMF models using simulated gene expression generated data. The simulated profiles were based on transcription factor (TF)-gene regulatory relationships inferred from breast tissue from GTEx using the Passing Attributes between Networks for Data Assimilation (PANDA) algorithm (Glass et al. 2013; Lonsdale et al. 2013). Gene expression profiles were then simulated using the inferred gene regulatory network as input to the Stochastic Gene Network Simulator (SGNSim) (Ribeiro and Lloyd-Price 2007; Tripathi et al. 2017). In addition, since the aim of GRACKLE is to identify subgroup-specific gene regulatory patterns, we isolated 5 network modules and systematically upregulated genes in these 5 modules in a set of samples to simulate subgroup-specific gene regulatory activation. 

We assessed the ability of GRACKLE and the other algorithms to decompose the gene expression data into latent variables that matched the simulated subgroups (i.e., gene module differential gene regulation) (Fig. 2). We randomly stratified the simulated profiles into a 70/30 training/testing split, then calculated the decomposed W and H matrices from the training data. We next projected the testing gene expression data into W using nonnegative least squares based on the trained H matrix. We calculated the accuracy of alignment to the subgroups based on whether the corresponding latent variable in the decomposed sample and gene matrices have the highest loadings for samples and genes of an assigned subgroup. For example, we consider a result to be accurate if the samples in subgroup A and the corresponding upregulated genes in module A have the highest loadings in LV1 of the decomposed sample and gene matrices W and H. We reported the average accuracy over all subgroups. The overall process is outlined in the schematic below: 

<img width="500" alt="image" src="https://github.com/user-attachments/assets/19529211-503c-4dfc-aab9-d944802b3b07" />


In the simulation studies, we evaluated performance over λ_1 (penalization for sample similarity, S_S) and λ_2 (penalization for gene similarity, S_G) values [0, 1] at an interval of 0.1. We tested GRACKLE over varying levels of background gene expression noise, decreased network modularity, and increased network transitivity (i.e., the percent of graph nodes involved in triangles). We benchmarked GRACKLE against three comparable algorithms: NMF, GNMF, and a prior informed graph regularized NMF inspired by netNMF-sc, which we call pr-GNMF. The NMF model served as a baseline and corresponds to λ_1 =0 and λ_2 =0. For GNMF, the affinity matrix for the gene regularization was calculated with k-nearest neighbors and the λ parameter, which affects the degree of graph regularization, was optimized using the parameters defined by Cai et al. (number of nearest neighbors = 5, λ = [1,10, 10^2, 10^3, 10^4] (Cai et al. 2011) ). For pr-GNMF, the same GRN used for GRACKLE was used for regularization of the gene similarity matrix. For both GNMF and pr-GNMF, graph regularization was performed for the maximum value of λ_2 tested. To avoid overfitting, we calculated the average performance over 100 iterations. 

