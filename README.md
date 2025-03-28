# GRACKLE

Here, we present a novel NMF approach that applies Graph Regulari-zation Across Contextual KnowLedgE (GRACKLE) to learn latent representations. GRACKLE constrains the unsupervised NMF model using both sample similarity and gene similarity matrices. The model is flexible to incorporate multiple forms of sample metadata, such as phenotype or subtype labels, as well as molecular relationships, such as gene-gene interactions or pathway annotation. 

A schematic of GRACKLE is shown below. 

<img width="400" alt="image" src="https://github.com/user-attachments/assets/bd3436bb-644b-486b-8f52-034d81ab54ff" />

Before matrix decomposition, the sample metadata is transformed into the sample similarity matrix $$S_S$$ by taking the dot product. Similarly, the gene interactions constitute the gene similarity matrix $$S_G$$. These similarity matrices are used to calculate the graph Laplacians $$L_S$$ and $$L_G$$, which are incorporated into the loss function as,

$$\min_{W,H} ||Y - WH||_F^2 + \lambda_1 \text{TR}(W^T L_S W) + \lambda_2 \text{TR}(H^T L_G H) \quad $$

where $$\lambda_1$$ and $$\lambda_2$$ are tuning parameters for graph regularization. Accordingly, the multiplicative update functions from Equation (4) are changed to,

$$w_{ik} \leftarrow w_{ik} \times \frac{YH^T + \lambda_1 S_S W}{WHH^T + \lambda_1 D_S W}, \quad h_{jk} \leftarrow h_{jk} \times \frac{W^T Y + \lambda_2 S_G H}{HWW^T + \lambda_2 D_G H} \quad $$

The model is iteratively trained until a stopping criterion of relative change in $$H$$ is less than $$1 \times 10^{-4}$$ between iterations or a maximum of 100 iterations. These model parameters may be fine tuned based on your scientific question.

## Prerequisites

*   **R:**  A statistical computing language and environment.
    *   Download and install the latest version of R from {Link: CRAN https://cran.r-project.org/}: 
        *   `https://cran.r-project.org/`

### Installing Required R Packages

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


