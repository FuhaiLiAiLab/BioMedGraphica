# BioMedGraphica

![BMG-logo](./Figures/BMG-logo.png)

<p align="center">
  <a href="https://github.com/FuhaiLiAiLab/BioMedGraphica">
    <img src="https://img.shields.io/badge/GitHub-BioMedGraphica-181717?logo=github" alt="GitHub Repo">
  </a>
  <a href="https://github.com/FuhaiLiAiLab/BioMedGraphica">
    <img src="https://img.shields.io/badge/WebUI-BioMedGraphica-red" alt="WebUI">
  </a>
  <a href="https://github.com/FuhaiLiAiLab/BioMedGraphica">
    <img src="https://img.shields.io/badge/Knowledge%20Graph-BioMedGraphica-blue" alt="Knowledge Graph">
  </a>
</p>
<p align="center">
  <a href="https://huggingface.co/datasets/FuhaiLiAiLab/BioMedGraphica">
    <img src="https://img.shields.io/badge/Hugging%20Face-Dataset-FFD21E?logo=huggingface" alt="Dataset">
  </a>
  <a href="https://www.biorxiv.org/content/10.1101/2024.12.05.627020v2">
    <img src="https://img.shields.io/badge/bioRxiv-Paper-6f42c1" alt="Paper">
  </a>
</p>

Artificial intelligence (AI) is transforming scientific discovery by leveraging its scalable capabilities to integrate and analyze large-scale datasets for knowledge mining. Foundation models—such as large language models (LLMs) and large vision models (LVMs)—are enabling general-purpose AI, yet they rely on highly structured data that is rarely available in the biomedical domain.

In contrast, biomedical data remains **fragmented**, with knowledge scattered across publications and inconsistent databases, each using diverse nomenclature systems. These discrepancies, spanning from genes to clinical traits, present significant challenges for **AI in Precision Health and Medicine (AI4PHM)**.

To address this, we introduce **BioMedGraphica**, an **all-in-one platform** for:
- Biomedical data integration across 43 databases
- Unified **Text-Attributed Knowledge Graph (TAKG)** generation
- Multi-omics and prior knowledge-driven graph AI applications

---

![Figure 1: BioMedGraphica Workflow](./Figures/Figure1.png)

---

## 📦 What is BioMedGraphica?

- **2,306,921 Entities**  
- **27,232,091 Relations**  
- **11 Entity Types** and **30 Harmonized Relation Types**  
- Fully attributed with unique IDs and **textual descriptions**
- Enables **zero-shot / few-shot knowledge discovery** via relation prediction
- Generates **graph AI-ready subgraphs** tailored to custom datasets

![BMG-statistics](./Figures/Figure0-1.png)

---

## 🧠 Use Cases

- Graph-based AI model training (e.g., GNNs, LLMs with graphs)
- Discovery of novel disease mechanisms and pathways
- Target and drug prioritization in AI4PHM
- Multi-omics signaling graph construction

---

## 🖥️ GUI & Software

BioMedGraphica offers an intuitive **Windows-based GUI**, allowing users to:
- Input multi-omics and clinical files
- Perform entity recognition (hard/soft matching)
- Construct user-specific knowledge-omic signaling graphs
- Export graph-ready `.npy` files for downstream AI modeling

### 🎥 Demo Video

https://github.com/user-attachments/assets/4c042d78-6be7-4596-83a6-ddf1439b5190

---

## 🗃️ Dataset Downloads

We recommend using the **[BioMedGraphica-Conn](https://huggingface.co/datasets/FuhaiLiAiLab/BioMedGraphica/tree/main/BioMedGraphica-Conn)** dataset, which excludes isolated nodes to support efficient graph training.

---

## 📁 Repository Structure

### 1. Data Collection

Check the structured data resource documentation [here](./DataCollection.md)

### 2. Data Processing

Scripts and processing logic are in the `BioMedGraphica-Raw` folder

### 3. GUI Software

The web app is now available online. Source code is available in [BioMedGraphica_WebApp](https://github.com/CallOfDady/BioMedGraphica_WebApp).

---

## 📚 Citation

If you use BioMedGraphica, please cite:

```bibtex
@article{zhang2024biomedgraphica,
  title={BioMedGraphica: An All-in-One Platform for Biomedical Prior Knowledge and Omic Signaling Graph Generation},
  author={Zhang, Heming and Liang, Shunning and Xu, Tim and Li, Wenyu and Huang, Di and Dong, Yuhan and Li, Guangfu and Miller, J Philip and Goedegebuure, S Peter and Sardiello, Marco and others},
  journal={bioRxiv},
  year={2024}
}
