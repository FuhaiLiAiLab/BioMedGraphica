# BioMedGraphica

## 1. How to Access and Reconstruct BioMedGraphica
**Data Access**: The datasets used in this project are available through APIs and links provided in the table. Ensure you have proper permissions and follow the terms of use for each source.

**Pre-Processing**: Use the provided Jupyter notebooks `.ipynb` to pre-process the data. Each notebook is designed to handle data retrieval, cleaning, and merging.
* RefSeq Raw Data Pre-Process
```python
# filter human data
def filter_data(input_file, output_file):

    df = pd.read_csv(input_file)
    
    filtered_df = df[(df['NCBI_tax_id'] == 9606) & (df['UniProtKB_tax_id'] == 9606)]
    
    filtered_df.to_csv(output_file, index=False)

# replace the input_file and output_file with the path of the files in your system
filter_data('gene_refseq_uniprotkb_collab', 'refseq_uniprot_human.csv')
```
**Output**: Processed data files `.csv` will be generated, which can be integrated into the BioMedGraphica framework.

## 2. API and Tool Integration
For entities that require API access,  the corresponding code is provided below and is also included in the respective `.ipynb` files to simplify data retrieval and integration.
* Ensembl API
```python
import pandas as pd
from pybiomart import Server

# List all available attributes
def list_attributes():
    server = Server(host='http://www.ensembl.org')
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    attributes = dataset.list_attributes()
    return attributes

attributes = list_attributes()

def fetch_ensembl_data(attributes):
    server = Server(host='http://www.ensembl.org')
    #https://www.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL
    #this link shows that hsapiens_gene_ensembl is the GRCh38.p14
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    
    response = dataset.query(attributes)
    
    return response

attributes = [] # choose the attributes you need
ensembl = fetch_ensembl_data(attributes)
```
* UniProt API
```python
import requests
from io import StringIO

def fetch_uniprot_data(params):
    url = "https://rest.uniprot.org/uniprotkb/stream"

    response = requests.get(url, params=params)

    if response.ok:
        tsv_data = StringIO(response.text)
        df = pd.read_csv(tsv_data, sep='\t')
        return df
    else:
        print("Failed to fetch data:", response.status_code)
        print(response.text)
        return None

params = {
        'fields': '', # choose the parameters you need
        'format': 'tsv',
        'query': '(model_organism:9606) AND (reviewed:true)',
        'sort': 'organism_name asc'
    }

uniprot = fetch_uniprot_data(params)
```
* CAS Data Retrieval
```python
import requests
import json

def fetch_and_write_annotations(base_url, total_pages, file_path):
    with open(file_path, 'w') as file:
        file.write('[')
        first_entry = True 
        
        for page in range(1, total_pages + 1):
            url = f"{base_url}?page={page}"
            print(f"Fetching data from: {url}") 
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                annotations = data.get('Annotations', {}).get('Annotation', [])
                
                for annotation in annotations:
                    if not first_entry:
                        file.write(',')
                    json.dump(annotation, file)
                    first_entry = False
            else:
                print(f"Failed to retrieve data for page {page}: {response.status_code}")
                continue
        
        file.write(']')

base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/CAS/JSON"
total_pages = 2789 # updata the total number if needed
file_path = "CAS_data.json"

fetch_and_write_annotations(base_url, total_pages, file_path)
```
* UNII Data Retrieval
```python
import requests
import json

def fetch_and_write_annotations(base_url, total_pages, file_path):
    with open(file_path, 'w') as file:
        file.write('[')
        first_entry = True 
        
        for page in range(1, total_pages + 1):
            url = f"{base_url}?page={page}"
            print(f"Fetching data from: {url}") 
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                annotations = data.get('Annotations', {}).get('Annotation', [])
                
                for annotation in annotations:
                    if not first_entry:
                        file.write(',')
                    json.dump(annotation, file)
                    first_entry = False
            else:
                print(f"Failed to retrieve data for page {page}: {response.status_code}")
                continue
        
        file.write(']')

base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/UNII/JSON"
total_pages = 153 # updata the total number if needed
file_path = "unii_data.json"

fetch_and_write_annotations(base_url, total_pages, file_path)
```
* KEGG Data Retrieval
```python
from bioservices import KEGG
import pandas as pd

k = KEGG()
k.organism = "hsa"

# human
pathway_ids = k.pathwayIds

def format_value(value):
    if isinstance(value, list):
        return ';'.join(format_value(item) for item in value)
    elif isinstance(value, dict):
        return ';'.join(f"{k}: {v}" for k, v in value.items())
    else:
        return str(value)

def format_for_dataframe(data):
    return {key: format_value(value) for key, value in data.items()}

full_pathway = pd.DataFrame()

for pid in pathway_ids:
    try:
        pathway_info = k.get(pid)  
        dict_data = k.parse(pathway_info)
        formatted_dict_data = format_for_dataframe(dict_data) 
        df = pd.DataFrame([formatted_dict_data])
        full_pathway = pd.concat([full_pathway, df], ignore_index=True)
    except Exception as e:
        print(f"Error processing pathway {pid}: {e}")

print(full_pathway.head())

full_pathway.to_csv('full_kegg_pathways.csv', index=False)
```
```R
# R code
library(graphite)
library(org.Hs.eg.db)
kps <- pathways("hsapiens", "kegg")
names_pathways <- names(kps)
n_p <- length(kps)
KeggPathways <- list()
for (i in 1:n_p){
  pt <- kps[[i]]
  pt <- convertIdentifiers(pt, 'symbol')
  et1 <- pt@mixedEdges
  et2 <- pt@protEdges
  et <- rbind(et1, et2)
  KeggPathways[[i]] <- et
}
names(KeggPathways) <- names_pathways
#
KeggGenes <- c()
for (i in 1:length(KeggPathways)){
  pt <- KeggPathways[[i]]
  gs <- c()
  xt <- pt[,1]; idxt <- regexpr('SYMBOL', xt); nt1 <- sum(idxt>0)
  if (nt1>0){gs1 <- pt[,2]; gs1 <- gs1[which(idxt>0)]; gs <- union(gs, gs1)}
  xt <- pt[,3]; idxt <- regexpr('SYMBOL', xt); nt2 <- sum(idxt>0)
  if (nt2>0){gs2 <- pt[,4]; gs2 <- gs2[which(idxt>0)]; gs <- union(gs, gs1)}
  if (max(nt1, nt2)>0){ KeggGenes <- union(KeggGenes, gs)}
}

# Step 3: Create KEGG edges list and add pathway IDs
KeggEdgesList <- c("source_type", "source", "target_type", "target", "direction", "edge_type", "pathway_name", "pathway_id")

for (i in 1:length(KeggPathways)){
  pathway_name <- names_pathways[[i]]
  pathway_id <- kps[[i]]@id # Extract pathway ID
  edgeTable <- KeggPathways[[i]]
  edgeTable <- cbind(edgeTable, rep(pathway_name, nrow(edgeTable)), rep(pathway_id, nrow(edgeTable))) # Add pathway name and ID
  print(dim(edgeTable))
  KeggEdgesList <- rbind(KeggEdgesList, edgeTable)
}

# Step 4: Finalize and save the list
colnames(KeggEdgesList) <- c("source_type", "source", "target_type", "target", "direction", "edge_type", "pathway_name", "pathway_id")
KeggEdgesList <- KeggEdgesList[-1,] # Remove the first row
unique_k <- unique(KeggEdgesList)
write.csv(unique_k, "full_kegg_pathway_list_with_id.csv")
```
* WikiPathway Data Retrieval
```python
import requests
import pandas as pd

url = "https://webservice.wikipathways.org/listPathways"
params = {
    "format": "json"
}

response = requests.get(url, params=params)
data = response.json()

pathways = data['pathways']

human_pathways = []

for pathway in pathways:
    if pathway['species'] == "Homo sapiens":
        human_pathways.append({
            'id': pathway['id'],
            'name': pathway['name'],
            'url': pathway['url'],
            'revision': pathway['revision']
        })

df = pd.DataFrame(human_pathways)
df.to_csv('human_pathways.csv', index=False)
```
* ChemIDplus Data Retrieval
```python
import pandas as pd
import json
import requests

def fetch_and_write_annotations(base_url, total_pages, file_path):
    with open(file_path, 'w') as file:
        file.write('[')
        first_entry = True 
        
        for page in range(1, total_pages + 1):
            url = f"{base_url}?page={page}"
            print(f"Fetching data from: {url}") 
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                annotations = data.get('Annotations', {}).get('Annotation', [])
                
                for annotation in annotations:
                    if not first_entry:
                        file.write(',')
                    json.dump(annotation, file)
                    first_entry = False
            else:
                print(f"Failed to retrieve data for page {page}: {response.status_code}")
                continue
        
        file.write(']')

base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/CAS/JSON"
total_pages = 2789  # Update to reflect the total number of pages
file_path = "CAS_data.json"

fetch_and_write_annotations(base_url, total_pages, file_path)

print(f"Data saved to {file_path}")
```
## 3. Detailed Steps of Data Download and Pre-processing of Entity
This section outlines the detailed steps for downloading and pre-processing data for each database used in Entity part. The procedures are designed to ensure data consistency and compatibility for downstream integration.
### 3.1 General Workflow
- **Data Source Identification**: Identify the official database or API for the required data.
- **Data Download**: Use appropriate tools or scripts to download the data.
- **Pre-processing**: Use custom scripts tailored to each database to clean, normalize, and restructure the data for integration.

### 3.2 Entity
Below are the key data sources used, along with their processing scripts and output files:
| Database | Download Link | Processing script | Output file |
| :-------: | :-------: | :-------: | :-------: |
| Ensembl | BioMart API | biomedgraphica_gene.ipynb | biomedgraphica_gene.csv |
|  | BioMart API | biomedgraphica_transcript.ipynb | biomedgraphica_transcript.csv |
|  | BioMart API | biomedgraphica_protein.ipynb | biomedgraphica_protein.csv |
| OMIM | [Link](https://omim.org/static/omim/data/mim2gene.txt) | biomedgraphica_gene.ipynb | biomedgraphica_gene.csv |
| HGNC  | [Link](https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit) | biomedgraphica_gene.ipynb | biomedgraphica_gene.csv |
| NCBI Gene | [Link1](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz) [Link2](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz) | biomedgraphica_gene.ipynb | biomedgraphica_gene.csv |
| RefSeq | [Link](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz) | biomedgraphica_gene.ipynb | biomedgraphica_gene.csv |
|  | [Link](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.3.summary.txt.gz) | biomedgraphica_transcript.ipynb | biomedgraphica_transcript.csv |
|  | [Link1](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz) [Link2](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz) | biomedgraphica_protein.ipynb | biomedgraphica_protein.csv |
| RNACentral | [Link](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/ensembl.tsv) | biomedgraphica_transcript.ipynb | biomedgraphica_transcript.csv |
| UniProt | API | biomedgraphica_protein.ipynb | biomedgraphica_protein.csv |
| Reactome | [Link](https://reactome.org/download/current/ReactomePathways.txt) | biomedgraphica_pathway.ipynb | biomedgraphica_pathway.csv |
| KEGG | Fetching data via R and Python | biomedgraphica_pathway.ipynb | biomedgraphica_pathway.csv |
| WikiPathways | Fetching data via Python | biomedgraphica_pathway.ipynb | biomedgraphica_pathway.csv |
| Pathway Ontology | [Link](https://download.rgd.mcw.edu/ontology/pathway/pathway.obo) | biomedgraphica_pathway.ipynb | biomedgraphica_pathway.csv |
| ComPath | [Link](https://compath.scai.fraunhofer.de/export_mappings) | biomedgraphica_pathway.ipynb | biomedgraphica_pathway.csv |
| HMDB | [Link](https://hmdb.ca/downloads) | biomedgraphica_metabolome.ipynb | biomedgraphica_metabolome.csv |
| ChEBI | [Link](https://www.ebi.ac.uk/chebi/chebiOntology.do?chebiId=77746) | biomedgraphica_metabolome.ipynb | biomedgraphica_metabolome.csv |
|  | [Link1](https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi_3star.tsv) [Link2](https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession_3star.tsv ) | biomedgraphica_drug.ipynb | biomedgraphica_drug.csv |
| NCBI Taxonomy | [Link](https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip) | biomedgraphica_microbiome.ipynb | biomedgraphica_microbiome.csv |
| SILVA | [Link1](https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/taxmap_embl-ebi_ena_lsu_ref_138.2.txt.gz) [Link2](https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/taxmap_embl-ebi_ena_ssu_ref_138.2.txt.gz) | biomedgraphica_microbiome.ipynb | biomedgraphica_microbiome.csv |
| Greengenes | [Link](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/greengenes.tsv) | biomedgraphica_microbiome.ipynb | biomedgraphica_microbiome.csv |
| RDP | [Link](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/rdp.tsv) | biomedgraphica_microbiome.ipynb | biomedgraphica_microbiome.csv |
| GTDB | [Link1](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz) [Link2](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz) | biomedgraphica_microbiome.ipynb | biomedgraphica_microbiome.csv |
| HPO | [Link](https://hpo.jax.org/data/ontology) | biomedgraphica_phenotype.ipynb | biomedgraphica_phenotype.csv |
| UMLS | [Link](https://download.nlm.nih.gov/umls/kss/2024AA/umls-2024AA-full.zip?_gl=1*14ig82q*_ga*MTA5NTI1Nzc2My4xNzEwOTU5NjM5*_ga_7147EPK006*MTcyMzU3NDM0NC41My4xLjE3MjM1NzUyNzYuMC4wLjA.*_ga_P1FPTH9PL4*MTcyMzU3NDM0NC41My4xLjE3MjM1NzUyNzYuMC4wLjA) | biomedgraphica_phenotype.ipynb | biomedgraphica_phenotype.csv |
|  | [Link](https://download.nlm.nih.gov/umls/kss/2024AA/umls-2024AA-full.zip?_gl=1*14ig82q*_ga*MTA5NTI1Nzc2My4xNzEwOTU5NjM5*_ga_7147EPK006*MTcyMzU3NDM0NC41My4xLjE3MjM1NzUyNzYuMC4wLjA.*_ga_P1FPTH9PL4*MTcyMzU3NDM0NC41My4xLjE3MjM1NzUyNzYuMC4wLjA) | biomedgraphica_disease.ipynb | biomedgraphica_disease.csv |
| ICD10 | [Link](https://icdcdn.who.int/static/releasefiles/2024-01/mapping.zip) | biomedgraphica_disease.ipynb | biomedgraphica_disease.csv |
| ICD11 | [Link](https://icdcdn.who.int/static/releasefiles/2024-01/SimpleTabulation-ICD-11-MMS-en.zip) | biomedgraphica_disease.ipynb | biomedgraphica_disease.csv |
| Disease Ontology | [Link](https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/DOreports/allXREFinDO.tsv) | biomedgraphica_disease.ipynb | biomedgraphica_disease.csv |
| MeSH | [Link](https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2024.xml) | biomedgraphica_disease.ipynb | biomedgraphica_disease.csv |
| SNOMED-CT | [Link](https://download.nlm.nih.gov/umls/kss/IHTSDO2024/IHTSDO20240801/SnomedCT_InternationalRF2_PRODUCTION_20240801T120000Z.zip?_gl=1*xret7k*_ga*MTA5NTI1Nzc2My4xNzEwOTU5NjM5*_ga_7147EPK006*MTcyMzU4ODA3OC41NC4xLjE3MjM1ODgyNDYuMC4wLjA.*_ga_P1FPTH9PL4*MTcyMzU4ODA3OS41NC4xLjE3MjM1ODgyNDYuMC4wLjA) | biomedgraphica_disease.ipynb | biomedgraphica_disease.csv |
| Mondo | [Link1](https://github.com/monarch-initiative/mondo/blob/master/reports/xrefs.tsv) [Link2](https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo) | biomedgraphica_disease.ipynb | biomedgraphica_disease.csv |
| PubChem | [Link1](https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22compound%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_compound_cache_n9g4FuJwh8yw5oX_B4fM1sZvYQ9Fdps_4RqAc_oLknL6Eq4%22,%22where%22:{%22ands%22:[{%22input%22:{%22type%22:%22netcachekey%22,%22idtype%22:%22cid%22,%22key%22:%22n9g4FuJwh8yw5oX_B4fM1sZvYQ9Fdps_4RqAc_oLknL6Eq4%22}}]}}) [Link2](https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22compound%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_compound_cache_66xMYpgP_bPKmf-Affi2qbwQOHAN461911i2McxJpDDMUJg%22,%22where%22:{%22ands%22:[{%22input%22:{%22type%22:%22netcachekey%22,%22idtype%22:%22cid%22,%22key%22:%2266xMYpgP_bPKmf-Affi2qbwQOHAN461911i2McxJpDDMUJg%22}}]}}) | biomedgraphica_drug.ipynb | biomedgraphica_drug.csv |
| CAS | Fetching data via Python | biomedgraphica_drug.ipynb | biomedgraphica_drug.csv |
| NDC | [Link](https://www.accessdata.fda.gov/cder/ndctext.zip) | biomedgraphica_drug.ipynb | biomedgraphica_drug.csv |
| UNII | [Link](https://precision.fda.gov/uniisearch/archive/latest/UNII_Data.zip); Fetching data via Python| biomedgraphica_drug.ipynb | biomedgraphica_drug.csv |
| DrugBank | [Link](https://go.drugbank.com/releases/5-1-12/downloads/all-drug-links) | biomedgraphica_drug.ipynb | biomedgraphica_drug.csv |
| CTD | [Link](https://ctdbase.org/reports/CTD_chemicals.csv.gz) | biomedgraphica_exposure.ipynb | biomedgraphica_exposure.csv |
| ToxCast | [Link](https://clowder.edap-cluster.com/files/6114f600e4b0856fdc65865c) | biomedgraphica_exposure.ipynb | biomedgraphica_exposure.csv |
| ChemIDplus | [Link](https://go.drugbank.com/releases/5-1-12/downloads/all-drug-links) | biomedgraphica_exposure.ipynb | biomedgraphica_exposure.csv |


### 3.3 Database-Specific Steps
#### 3.3.1 Ensembl

Ensembl data is used across three entities: `gene`, `transcript`, and `protein`. Data is downloaded using the API described in the `API and Tool Integration` section, and specific attributes are retained for each entity during pre-processing. Custom processing scripts are utilized to ensure data consistency and relevance.

##### **Data Download**
- Use the [BioMart](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-10-22) API to fetch data for `gene`, `transcript`, and `protein` entities.
- The API allows querying specific datasets and attributes.

##### **Pre-processing**
- Custom Ensembl section in each scripts (`biomedgraphica_gene.ipynb`, `biomedgraphica_transcript.ipynb`, `biomedgraphica_protein.ipynb`) are used to clean and format the data for integration.
- For `gene` entity, retain the following attributes:  
```python
attributes=['ensembl_gene_id', 'ensembl_gene_id_version','start_position', 'end_position','gene_biotype', 'hgnc_id', 'hgnc_symbol']
```
- For `transcript` entity, retain the following attributes:
```python
attributes=['ensembl_transcript_id' ,'ensembl_transcript_id_version', 'ensembl_gene_id', 'external_gene_name', 'external_transcript_name', 'transcript_biotype', 'refseq_mrna', 'refseq_ncrna', 'transcript_mane_select']
```
- For `protein` entity, retain the following attributes:
```python
attributes=['ensembl_peptide_id', 'ensembl_peptide_id_version','uniprotswissprot', 'refseq_peptide', 'entrezgene_id', 'external_gene_name']
```
---
#### 3.3.2 OMIM - Need Registration
The OMIM database provides two key data sources used in this project: `mim2gene.txt` and `genemap2.txt`. These datasets are essential for gene-level integration and are pre-processed using the `biomedgraphica_gene.ipynb`.
##### **Data Download**
- **mim2gene.txt**:
   - This file can be directly downloaded without registration.
   - It contains mappings between OMIM IDs and gene IDs.

- **genemap2.txt**:
   - This file requires registration on the [OMIM downloads page](https://omim.org/downloads).
   - Note: This dataset is subject to non-commercial use restrictions. Please review the OMIM terms of use carefully before downloading and processing.
##### **Pre-processing**
- The pre-processing for both files is implemented in the `biomedgraphica_gene.ipynb` notebook.
---
#### 3.3.3 HGNC

The HGNC database provides essential gene-level information that is used in this project. Data is retrieved using HGNC's Biomart tool, which allows precise filtering and extraction of relevant fields. The processed data is then used for gene-level integration.
##### **Data Download**
- Access the data via the `Download link` in the `Entity` section below.
- **Fields Selected**:
  - `HGNC ID`
  - `Approved symbol`
  - `Approved name`
  - `NCBI Gene ID`
  - `Ensembl gene ID`
##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_gene.ipynb` notebook.
---
#### 3.3.4 NCBI Gene
The NCBI Gene database provides critical gene-level information, particularly for human data and mapping relationships between NCBI Gene IDs and Ensembl IDs. Two files, `gene_info` and `gene2ensembl`, are utilized for this purpose.

##### **Data Download**
- The files can be accessed via the `Download link` in the `Entity` section of this document.
- **gene_info**:
   - Contains comprehensive information on genes, including descriptions, synonyms, and taxonomy.
   - Only human data is extracted for this project.

- **gene2ensembl**:
   - Provides a mapping between NCBI Gene IDs and Ensembl IDs.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_gene.ipynb` notebook.
---
#### 3.3.5 RefSeq

The RefSeq database is utilized for `gene`, `transcript`, and `protein` entities in this project. Different files are used for each entity, and pre-processing is performed in the corresponding Jupyter Notebooks.
##### **Data Download**
- The following files are used for each entity and can be accessed via the `Download link` in the `Entity` section of this document:
- **Gene**:
   - File: `gene_refseq_uniprotkb_collab`
   - Due to the large size of the file, only human data is extracted using the **RefSeq Data Pre-Process code** in `biomedgraphica_gene.ipynb`.

- **Transcript**:
   - File: `MANE.GRCh38.v1.3.summary.txt`

- **Protein**:
   - Files: `refseq_uniprot_human.csv` and `gene2ensembl`
   - Human-specific data from `refseq_uniprot_human.csv` is extracted using the **RefSeq Data Pre-Process code** in `biomedgraphica_protein.ipynb`.

##### **Pre-processing**
- Pre-processing is performed using the following notebooks:
  - **Gene-level data**: `biomedgraphica_gene.ipynb`
    - Extracts human-specific data from `gene_refseq_uniprotkb_collab` and formats it for integration.
  - **Transcript-level data**: `biomedgraphica_gene.ipynb`
    - Directly processes `MANE.GRCh38.v1.3.summary.txt`.
  - **Protein-level data**: `biomedgraphica_protein.ipynb`
    - Extracts human-specific data from `refseq_uniprot_human.csv` and processes `gene2ensembl`.
---
#### 3.3.6 RNACentral

The RNACentral database provides comprehensive transcript data for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration or filtering is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_transcript.ipynb` notebook.
---
#### 3.3.7 UniProt

The UniProt database is used to retrieve protein-level data for this project. Data is downloaded programmatically using the UniProt API.

##### **Data Download**
- Data is accessed via the UniProt API.
- The specific API query used:
```python
params = {
        'fields': 'accession,protein_name,gene_primary,xref_ensembl_full,xref_geneid',
        'format': 'tsv',
        'query': '(model_organism:9606) AND (reviewed:true)',
        'sort': 'organism_name asc'
    }
```
- Example API usage and detailed download code can be found in the `API and Tool Integration` section or in the `biomedgraphica_protein.ipynb` notebook.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_protein.ipynb` notebook.
---
#### 3.3.8 Reactome

The Reactome database provides pathway-level data for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_pathway.ipynb` notebook.
---
#### 3.3.9 KEGG
The KEGG database provides pathway and gene-level data for this project. Data is obtained using both R and Python, depending on the specific use case.

##### **Data Download**
- KEGG data is retrieved programmatically using Python and R. Example API usage and detailed download code can be found in the `API and Tool Integration` section or in the `biomedgraphica_pathway.ipynb` notebook.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_pathway.ipynb` notebook.
---
#### 3.3.10 WikiPathways

The WikiPathways database provides pathway-level data for this project. Data is retrieved programmatically using the WikiPathways API.

##### **Data Download**
- Data is accessed via the WikiPathways API.
- Example API usage and detailed download code can be found in the `API and Tool Integration` section or in the `biomedgraphica_pathway.ipynb` notebook.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_pathway.ipynb` notebook.
---
#### 3.3.11 Pathway Ontology
The Pathway Ontology database provides hierarchical pathway data for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_pathway.ipynb` notebook.
---
#### 3.3.12 ComPath
The ComPath database provides ID mapping between KEGG, WikiPathways, and Pathway Ontology, facilitating pathway data integration.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_pathway.ipynb` notebook.
---
#### 3.3.13 HMDB

The HMDB database provides comprehensive metabolite information for this project.

##### **Data Download**
- Click the `Download link` in the `Entity` section of this document to download the "All Metabolites" XML file.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_metabolome.ipynb` notebook.
---
#### 3.3.14 ChEBI

The ChEBI database is used for both `metabolite` and `drug` entities in this project. It provides curated information about small molecular entities.

##### **Data Download**
- **Metabolite**:
   - The provided download link directs to a page with human metabolites.
   - Click the "Download as Tab-delimited" option on the page to download the file.

- **Drug**:
   - Click the `Download link` in the `Entity` section of this document to download the data.
   - No registration is required.

##### **Pre-processing**
- Pre-processing is performed using the following notebooks:
  - **Metabolite data**: `biomedgraphica_metabolome.ipynb`
  - **Drug data**: `biomedgraphica_drug.ipynb`
---
#### 3.3.15 NCBI Taxonomy

The NCBI Taxonomy database provides hierarchical and taxonomic data for the `microbiota` entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_microbiome.ipynb` notebook.
---
#### 3.3.16 SILVA

The SILVA database provides curated and comprehensive microbiota data for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_microbiome.ipynb` notebook.

---
#### 3.3.17 Greengenes
The Greengenes database offers taxonomic data for the microbiota entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_microbiome.ipynb` notebook.

---
#### 3.3.18 RDP

The RDP (Ribosomal Database Project) database provides curated microbiota-related data for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_microbiome.ipynb` notebook.

---

#### 3.3.19 GTDB

The GTDB (Genome Taxonomy Database) offers genomic and taxonomic data for the microbiota entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_microbiome.ipynb` notebook.
---
#### 3.3.20 HPO (Human Phenotype Ontology)

The HPO database provides standardized phenotype terminology for the `phenotype` entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document. Download the `HP.OBO` file.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_phenotype.ipynb` notebook.
- The notebook cleans up phenotype labels by:
  - Identifying and systematically removing predefined descriptive expressions.
  - Combining all duplicate entries to ensure data consistency.
---

#### 3.3.21 UMLS - Need Registration

The UMLS database provides comprehensive mappings of medical terms, including phenotype-related data, for this project.

##### **Data Download**
- Registration is required to access the UMLS data. Please follow the instructions provided on the UMLS website to complete the registration process and download the necessary files.

##### **Pre-processing**
- Pre-processing is performed using the following notebooks:
  - **Phenotype data**: `biomedgraphica_phenotype.ipynb`
  - **Disease data**: `biomedgraphica_disease.ipynb`
---

#### 3.3.22 ICD-10

The ICD-10 database provides standardized disease classifications for the `disease` entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_disease.ipynb` notebook.

---

#### 3.3.23 ICD-11

The ICD-11 database offers updated disease classifications for the `disease` entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_disease.ipynb` notebook.
---

#### 3.3.24 Disease Ontology

The Disease Ontology database provides hierarchical disease terms and mappings for the `disease` entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_disease.ipynb` notebook.
---

#### 3.3.25 MeSH

The MeSH database provides hierarchical medical terms and classifications used in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_disease.ipynb` notebook.
---

#### 3.3.26 SNOMED-CT - Need Registration

The SNOMED-CT database provides detailed clinical terminology used in this project.

##### **Data Download**
- Registration is required to access the SNOMED-CT data. Please follow the instructions provided on the SNOMED-CT or UMLS registration system to complete the registration process and download the necessary files.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_disease.ipynb` notebook.

---
#### 3.3.27 MONDO

The MONDO database provides integrated disease ontology terms for the `disease` entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_disease.ipynb` notebook.

---
#### 3.3.28 PubChem

The PubChem database provides comprehensive chemical information used in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_drug.ipynb` notebook.
---

#### 3.3.29 CAS

The CAS database provides unique identifiers and annotations for chemical substances used in this project.

##### **Data Download**
- Data is fetched programmatically using a Python script that utilizes the PubChem API.
- Example Python script for fetching CAS data:
  ```python
  import requests
  import json

  def fetch_and_write_annotations(base_url, total_pages, file_path):
      with open(file_path, 'w') as file:
          file.write('[')
          first_entry = True 
          
          for page in range(1, total_pages + 1):
              url = f"{base_url}?page={page}"
              print(f"Fetching data from: {url}") 
              response = requests.get(url)
              if response.status_code == 200:
                  data = response.json()
                  annotations = data.get('Annotations', {}).get('Annotation', [])
                  
                  for annotation in annotations:
                      if not first_entry:
                          file.write(',')
                      json.dump(annotation, file)
                      first_entry = False
              else:
                  print(f"Failed to retrieve data for page {page}: {response.status_code}")
                  continue
          
          file.write(']')

  base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/CAS/JSON"
  total_pages = 2789  # Update the total number if needed
  file_path = "CAS_data.json"

  fetch_and_write_annotations(base_url, total_pages, file_path)
  ```
##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_drug.ipynb` notebook.
---
#### 3.3.30 NDC (National Drug Code)

The NDC database provides standardized drug identification codes for the `drug` entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_drug.ipynb` notebook.
---

#### 3.3.31 UNII (Unique Ingredient Identifier)

The UNII database provides unique identifiers for drug ingredients, split into two parts.

##### **Data Download**
1. **Part 1**:
   - Can be directly downloaded via the `Download link` in the `Entity` section of this document.

2. **Part 2**:
   - Data is fetched programmatically using a Python script via the PubChem API.
   - Example Python script:
     ```python
     import requests
     import json

     def fetch_and_write_annotations(base_url, total_pages, file_path):
         with open(file_path, 'w') as file:
             file.write('[')
             first_entry = True 
             
             for page in range(1, total_pages + 1):
                 url = f"{base_url}?page={page}"
                 print(f"Fetching data from: {url}") 
                 response = requests.get(url)
                 if response.status_code == 200:
                     data = response.json()
                     annotations = data.get('Annotations', {}).get('Annotation', [])
                     
                     for annotation in annotations:
                         if not first_entry:
                             file.write(',')
                         json.dump(annotation, file)
                         first_entry = False
                 else:
                     print(f"Failed to retrieve data for page {page}: {response.status_code}")
                     continue
             
             file.write(']')

     base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/UNII/JSON"
     total_pages = 153  # Update the total number if needed
     file_path = "unii_data.json"

     fetch_and_write_annotations(base_url, total_pages, file_path)
     ```

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_drug.ipynb` notebook.

---

#### 3.3.32 DrugBank

The DrugBank database provides detailed drug information for the `drug` entity in this project.

##### **Data Download**
- Registration for an academic account is required to access the DrugBank data.
- Once registered, data can be downloaded from the DrugBank website.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_drug.ipynb` notebook.
---
#### 3.3.33 CTD

The CTD database provides curated data on chemical exposures, diseases, and gene interactions for the `exposure` entity in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_exposure.ipynb` notebook.
---

#### 3.3.34 ToxCast

The ToxCast database provides high-throughput screening data for chemical toxicity used in this project for the `exposure` entity.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_exposure.ipynb` notebook.
---
#### 3.3.35 ChemIDplus

The ChemIDplus database provides curated chemical information and is used in this project for the `exposure` and `drug` entities.

##### **Data Download**
- Data is fetched programmatically using a Python script via the PubChem API.
- The following script retrieves the data and saves it as `CAS_data.json`:
  ```python
  import pandas as pd
  import json
  import requests

  def fetch_and_write_annotations(base_url, total_pages, file_path):
      with open(file_path, 'w') as file:
          file.write('[')
          first_entry = True 
          
          for page in range(1, total_pages + 1):
              url = f"{base_url}?page={page}"
              print(f"Fetching data from: {url}") 
              response = requests.get(url)
              if response.status_code == 200:
                  data = response.json()
                  annotations = data.get('Annotations', {}).get('Annotation', [])
                  
                  for annotation in annotations:
                      if not first_entry:
                          file.write(',')
                      json.dump(annotation, file)
                      first_entry = False
              else:
                  print(f"Failed to retrieve data for page {page}: {response.status_code}")
                  continue
          
          file.write(']')

  base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/CAS/JSON"
  total_pages = 2789  # Update to reflect the total number of pages
  file_path = "CAS_data.json"

  fetch_and_write_annotations(base_url, total_pages, file_path)

  print(f"Data saved to {file_path}")
  ```
##### **Pre-processing**
- Pre-processing is performed using the `biomedgraphica_exposure.ipynb` notebook.
---


## 4. Detailed Steps of Data Download and Pre-processing of Relation
This section outlines the databases used for relation extraction and mapping, detailing the download and pre-processing steps for integrating them into the knowledge graph.

### 4.1 Relation
Below are the key data sources used, along with their processing scripts and output files:
| Database | Download Link | Processing script | Output file |
| :-------: | :-------: | :-------: | :-------: |
| - | - | promoter-gene.ipynb | biomedgraphica_promoter_gene.csv |
| Ensembl | BioMart API | gene-transcript.ipynb | biomedgraphica_gene_transcript.csv |
| RefSeq | [Link](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene) | gene-transcript.ipynb | biomedgraphica_gene_transcript.csv |
| Ensembl | BioMart API | transcript-protein.ipynb | biomedgraphica_transcript_protein.csv |
| UniProt | API | transcript-protein.ipynb | biomedgraphica_transcript_protein.csv |
| RefSeq | [Link](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene) | transcript-protein.ipynb | biomedgraphica_transcript_protein.csv |
| BioGrid | [Link](https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.237/BIOGRID-ALL-4.4.237.mitab.zip) | protein-protein.ipynb | biomedgraphica_protein_protein.csv |
| STRING | [Link](https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz) | protein-protein.ipynb | biomedgraphica_protein_protein.csv |
| KEGG | Fetching data via R | protein-protein.ipynb | biomedgraphica_protein_protein.csv |
| KEGG | Fetching data via R | protein-pathway.ipynb | biomedgraphica_protein_pathway.csv |
| HPO | [Link](https://hpo.jax.org/data/annotations) | protein-phenotype.ipynb | biomedgraphica_protein_phenotype.csv |
| UniProt | [Link](https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Ccc_disease&format=tsv&query=%28*%29+AND+%28model_organism%3A9606%29) | protein-disease.ipynb | biomedgraphica_protein_disease.csv |
| DISEASES | [Link](https://download.jensenlab.org/human_disease_benchmark.tsv) | protein-disease.ipynb | biomedgraphica_protein_disease.csv |
| HPO | [Link](https://hpo.jax.org/data/annotations) | protein-disease.ipynb | biomedgraphica_protein_disease.csv |
| DisGeNet | API | protein-disease.ipynb | biomedgraphica_protein_disease.csv |
| KEGG | Fetching data via R | pathway-protein.ipynb | biomedgraphica_pathway_protein.csv |
| CTD | [Link](https://ctdbase.org/reports/CTD_chem_pathways_enriched.csv.gz) | pathway-exposure.ipynb | biomedgraphica_pathway_exposure.csv |
| KEGG | Fetching data via R | pathway-drug.ipynb | biomedgraphica_pathway_drug.csv |
| HMDB | [Link](https://hmdb.ca/downloads) | metabolite-protein.ipynb | biomedgraphica_metabolome_protein.csv |
| MetaNetX | [Link1](https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv) [Link2](https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_isom.tsv) | metabolite-metabolite.ipynb | biomedgraphica_metabolome_metabolite.csv |
| HMDB | [Link](https://hmdb.ca/downloads) | metabolite-disease.ipynb | biomedgraphica_metabolome_disease.csv |
| DisBiome | [Link](https://disbiome.ugent.be/export) | microbiota-disease.ipynb | biomedgraphica_microbiome_disease.csv |
| MDAD | [Link](https://github.com/Sun-Yazhou/MDAD/blob/master/MDAD.zip) | microbiota-drug.ipynb | biomedgraphica_microbiome_drug.csv |
| PharmacoMicrobiomics | [Link](http://pharmacomicrobiomics.com/view/relation/) | microbiota-drug.ipynb | biomedgraphica_microbiome_drug.csv |
| CTD | [Link](https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz) | exposure-gene.ipynb | biomedgraphica_exposure_gene.csv |
| CTD | [Link](https://ctdbase.org/reports/CTD_chem_pathways_enriched.csv.gz) | exposure-pathway.ipynb | biomedgraphica_exposure_pathway.csv |
| CTD | [Link](https://ctdbase.org/reports/CTD_chemicals_diseases.csv.gz) | exposure-disease.ipynb | biomedgraphica_exposure_disease.csv |
| HPO | [Link](https://hpo.jax.org/data/ontology) | phenotype-phenotype.ipynb | biomedgraphica_phenotype_phenotype.csv |
| HPO | [Link](https://hpo.jax.org/data/annotations) | phenotype-disease.ipynb | biomedgraphica_phenotype_disease.csv |
| HPO | [Link](https://hpo.jax.org/data/annotations) | disease-phenotype.ipynb | biomedgraphica_disease_phenotype.csv |
| DO | [Link](https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/refs/heads/main/src/ontology/HumanDO.obo) | disease-disease.ipynb | biomedgraphica_disease_disease.csv |
| DrugBank | [Link](https://go.drugbank.com/releases/5-1-12/downloads/target-all-polypeptide-ids) | drug-protein.ipynb | biomedgraphica_drug_protein.csv |
| BindingDB | [Link](https://www.bindingdb.org/bind/downloads/BindingDB_All_202409_tsv.zip) | drug-protein.ipynb | biomedgraphica_drug_protein.csv |
| DrugCentral | [Link](https://drugcentral.org/ActiveDownload) | drug-protein.ipynb | biomedgraphica_drug_protein.csv |
| KEGG | Fetching data via R | drug-pathway.ipynb | biomedgraphica_drug_pathway.csv |
| HMDB | [Link](https://hmdb.ca/downloads) | drug-metabolite.ipynb | biomedgraphica_drug_metabolome.csv |
| MDAD | [Link](https://github.com/Sun-Yazhou/MDAD/blob/master/MDAD.zip) | drug-microbiota.ipynb | biomedgraphica_drug_microbiota.csv |
| PharmacoMicrobiomics | [Link](http://pharmacomicrobiomics.com/view/relation/) | drug-microbiota.ipynb | biomedgraphica_drug_microbiota.csv |
| SIDER | [Link](http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz) | drug-phenotype.ipynb | biomedgraphica_drug_phenotype.csv |
| DrugCentral | [Link](https://drugcentral.org/ActiveDownload) | drug-disease.ipynb | biomedgraphica_drug_disease.csv |
| DrugBank | [Link](https://go.drugbank.com/releases/5-1-12/downloads/all-full-database) | drug-drug.ipynb | biomedgraphica_drug_drug.csv |

### 4.2 Database-Specific Steps
#### 4.2.1 Ensembl

The Ensembl database provides data for gene-transcript and transcript-protein relations, retrieved using the [BioMart](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-10-22).

##### **Data Download**
- Data is accessed via the [BioMart](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-10-22), specifying the attributes relevant to each relation type:
  - **Gene-Transcript Relation**:
     - Attributes: `['ensembl_gene_id', 'ensembl_transcript_id']`

  - **Transcript-Protein Relation**:
     - Attributes: `['ensembl_transcript_id', 'ensembl_protein_id']`

##### **Pre-processing**
- **Gene-Transcript Relation**:
   - Pre-processing is performed using the `gene-transcript.ipynb` notebook.
- **Transcript-Protein Relation**:
   - Pre-processing is performed using the `transcript-protein.ipynb` notebook.
---

#### 4.2.2 RefSeq

The RefSeq database provides data for `gene-transcript` and `transcript-protein` relations in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- **Gene-Transcript Relation**:
   - Pre-processing is performed using the `gene-transcript.ipynb` notebook.

- **Transcript-Protein Relation**:
   - Pre-processing is performed using the `transcript-protein.ipynb` notebook.
---

#### 4.2.3 UniProt

The UniProt database provides data for `transcript-protein` and `protein-disease` relations in this project.

##### **Data Download**
- Data is retrieved programmatically using the UniProt API with the following parameters:
  ```python
  params = {
      'fields': 'accession,xref_ensembl',  # choose the parameters you need
      'format': 'tsv',
      'query': '(model_organism:9606) AND (reviewed:true)',  # human reviewed entries
      'sort': 'organism_name asc'
  }
  ```
##### **Pre-processing**
- **Transcript-Protein Relation:**:
   - Pre-processing is performed using the `transcript-protein.ipynb` notebook.

- **Protein-Disease Relation**:
   - Pre-processing is performed using the `protein-disease.ipynb` notebook.
---

#### 4.2.4 BioGRID

The BioGRID database provides curated protein-protein interaction data for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `protein-protein.ipynb` notebook.
---

#### 4.2.5 STRING

The STRING database provides protein-protein interaction data, including interaction confidence scores.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `protein-protein.ipynb` notebook.
---

#### 4.2.6 KEGG

The KEGG database provides data for multiple relations in this project, including `protein-protein`, `protein-pathway`, `pathway-protein`, `pathway-drug`, and `drug-pathway`.

##### **Data Download**
- KEGG data is retrieved using R scripts. The specific code can be found in the `Entity` section of this document.

##### **Pre-processing**
- **Protein-Protein Relation**:
   - Pre-processing is performed using the `protein-protein.ipynb` notebook.

- **Protein-Pathway Relation**:
   - Pre-processing is performed using the `protein-pathway.ipynb` notebook.

- **Pathway-Protein Relation**:
   - Pre-processing is performed using the `pathway-protein.ipynb` notebook.

- **Pathway-Drug Relation**:
   - Pre-processing is performed using the `pathway-drug.ipynb` notebook.

- **Drug-Pathway Relation**:
   - Pre-processing is performed using the `drug-pathway.ipynb` notebook.
---

#### 4.2.7 HPO (Human Phenotype Ontology)

The HPO database provides standardized phenotype information and is used for multiple relations in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- **Protein-Phenotype Relation**:
   - Pre-processing is performed using the `protein-phenotype.ipynb` notebook.

- **Protein-Disease Relation**:
   - Pre-processing is performed using the `protein-disease.ipynb` notebook.

- **Phenotype-Phenotype Relation**:
   - Pre-processing is performed using the `phenotype-phenotype.ipynb` notebook.

- **Phenotype-Disease Relation**:
   - Pre-processing is performed using the `phenotype-disease.ipynb` notebook.

- **Disease-Phenotype Relation**:
   - Pre-processing is performed using the `disease-phenotype.ipynb` notebook.
---

#### 4.2.8 DISEASES
The DISEASES database provides curated associations between proteins and diseases for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `protein-disease.ipynb` notebook.
---
#### 4.2.9 DisGeNET

The DisGeNET database provides curated gene-disease associations for this project.

##### **Data Download**
- Data is retrieved programmatically via the DisGeNET API.
- Registration for an academic account is required to access the API. Follow the instructions on the DisGeNET website to create an account and obtain your API key.
- Example Python script for fetching gene-disease data:
  ```python
  import pandas as pd
  import requests
  import json
  import time
  import os
  from tqdm import tqdm  
  from requests.exceptions import ConnectionError  
  import warnings
  from urllib3.exceptions import InsecureRequestWarning

  warnings.filterwarnings('ignore', category=InsecureRequestWarning)

  API_KEYS = ["your_api_key_here"]  # Replace with your valid API key
  file_path = 'ncbi_gene_id_with_symbol.csv'  # Path to input file
  output_file = "gene_disease_umls.csv"  # Output file

  def send_request(gene_id, page_number, api_key):
      params = {"gene_ncbi_id": str(gene_id), "page_number": str(page_number)}
      headers = {"Authorization": api_key, "accept": "application/json"}
      response = requests.get("https://api.disgenet.com/api/v1/gda/summary", params=params, headers=headers, verify=False)
      return response.json() if response.ok else None

  # Example processing function
  def process_data():
      for chunk in pd.read_csv(file_path, chunksize=100):
          for _, row in chunk.iterrows():
              gene_id = row["GeneID"]
              page_number = 0
              while True:
                  data = send_request(gene_id, page_number, API_KEYS[0])
                  if data and "payload" in data:
                      for association in data["payload"]:
                          print(f"Gene: {gene_id}, Disease: {association.get('diseaseName')}")
                          # Write to output
                  page_number += 1
                  if page_number >= data.get("paging", {}).get("totalPages", 1):
                      break

  process_data()
  ```
##### **Pre-processing**
- Pre-processing is performed using the `gene-disease.ipynb` notebook.
---
#### 4.2.10 HMDB (Human Metabolome Database)

The HMDB database provides curated metabolite data and is used for multiple relations in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- **Metabolite-Protein Relation**:
   - Pre-processing is performed using the `metabolite-protein.ipynb` notebook.

- **Metabolite-Disease Relation**:
   - Pre-processing is performed using the `metabolite-disease.ipynb` notebook.

- **Drug-Metabolite Relation**:
   - Pre-processing is performed using the `drug-metabolite.ipynb` notebook.

- The notebooks extract the relevant relationships from the downloaded HMDB data, standardize the formats, and integrate them into the knowledge graph.
---

#### 4.2.11 DisBiome

The DisBiome database provides curated microbiota-disease associations for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `microbiota-disease.ipynb` notebook.
---
#### 4.2.12 DrugBank

The DrugBank database provides comprehensive drug-related information and is used for multiple relations in this project.

##### **Data Download**
- Registration for an academic account is required to access the DrugBank data.
- Follow the instructions on the DrugBank website to create an account and download the necessary files.

##### **Pre-processing**
- **Drug-Protein Relation**:
   - Pre-processing is performed using the `drug-protein.ipynb` notebook.

- **Drug-Drug Relation**:
   - Pre-processing is performed using the `drug-drug.ipynb` notebook.
---

#### 4.2.13 DrugCentral

The DrugCentral database provides curated drug-related information and is used for multiple relations in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- The downloaded file is an SQL database containing all DrugCentral data.

##### **Pre-processing**
- **Drug-Protein Relation**:
   - Pre-processing is performed using the `drug-protein.ipynb` notebook.
   - Extracting the necessary relationships (e.g., drug-protein interactions) using SQL queries.

- **Drug-Disease Relation**:
   - Pre-processing is performed using the `drug-disease.ipynb` notebook.
   - Extracting the necessary relationships (e.g., drug-disease associations) using SQL queries.
---

#### 4.2.14 BindingDB

The BindingDB database provides data on molecular interactions, including drug-target binding affinities, and is used in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `drug-protein.ipynb` notebook.
- The downloaded file contains a large dataset. A Python script is used to filter the data for human-related entries to reduce the size of the dataset:
   ```python
    import pandas as pd

    def filter_human_bindingdb(input_file, output_file, chunksize=10000):
        df_bindingdb_human = pd.DataFrame()

        for i, chunk in enumerate(pd.read_csv(input_file, sep='\t', on_bad_lines='skip', iterator=True, chunksize=chunksize)):
            chunk_human = chunk[chunk['Target Source Organism According to Curator or DataSource'] == 'Homo sapiens']

            chunk_human.to_csv(output_file, mode='a', header=(i == 0), sep='\t', index=False)

    if __name__ == "__main__":
        input_file = 'BindingDB_All.tsv'  # Modify with your input file path
        output_file = 'BindingDB_Human.tsv'  # Modify with your desired output file path
        filter_human_bindingdb(input_file, output_file)
   ```

#### 4.2.15 SIDER

The SIDER database provides curated drug-phenotype (side effect) associations and is used in this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `drug-phenotype.ipynb` notebook.
---

#### 4.2.16 MetaNetX

The MetaNetX database provides curated metabolite-metabolite associations for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `metabolite-metabolite.ipynb` notebook.
---

#### 4.2.17 MDAD

The MDAD database provides curated microbiota-drug interaction data for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- **Microbiota-Drug Relation**:
   - Pre-processing is performed using the `microbiota-drug.ipynb` notebook.

- **Drug-Microbiota Relation**:
   - Pre-processing is performed using the `drug-microbiota.ipynb` notebook.
---

#### 4.2.18 PharmacoMicrobiomics

The PharmacoMicrobiomics database provides curated data on microbiota-drug interactions for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- **Microbiota-Drug Relation**:
   - Pre-processing is performed using the `microbiota-drug.ipynb` notebook.

- **Drug-Microbiota Relation**:
   - Pre-processing is performed using the `drug-microbiota.ipynb` notebook.
---

#### 4.2.19 DO (Disease Ontology)

The Disease Ontology (DO) database provides curated disease-disease associations and hierarchical relationships for this project.

##### **Data Download**
- The data can be accessed via the `Download link` in the `Entity` section of this document.
- No registration is required for downloading.

##### **Pre-processing**
- Pre-processing is performed using the `disease-disease.ipynb` notebook.
