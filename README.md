# GRACKLE: Graph Regularized NMF Across Contextual Knowledge
**Motivation**: Dysregulation of established transcriptomic regulation drive diseases and chronic conditions. However, identifying disease-specific gene signatures can be challenging due to the presence of
multiple co-occurring conditions and limited sample sizes. Unsupervised representation learning methods, such as matrix decomposition and deep learning, reduce high-dimensional data to interpretable
patterns yet lack explicit biological explainability. Incorporating prior biological knowledge directly into models can enhance understanding and address small sample sizes. Nevertheless, current models do
not jointly consider prior knowledge of molecular interactions and sample labels. 

**Results**: We present GRACKLE, a novel non-negative matrix factorization approach that applies
Graph Regularization Across Contextual KnowLedgE. GRACKLE integrates sample similarity and
gene similarity matrices based on sample metadata and molecular relationships, respectively. Simulation
studies show GRACKLE outperformed other NMF algorithms, especially with increased background
noise. GRACKLE effectively stratified breast tumor samples and identified condition-enriched
subgroups in individuals with Down syndrome. The model's latent representations aligned with known
biological patterns, such as autoimmune conditions and sleep apnea in Down syndrome. GRACKLE's
flexibility allows application to various data modalities, offering a robust solution for identifying context-
specific molecular mechanisms in biomedical research.


<img src="https://github.com/user-attachments/assets/ae8be5ae-fc1f-4d48-912a-3277ff14d97d" width="675" height="300">

# Download
Download code and data to run analysis via:
```
git clone https://github.com/lagillenwater/GRACKLE
```
