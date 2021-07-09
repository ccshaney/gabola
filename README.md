# GABOLA: A Reliable Gap-Filling Strategy for de novo Chromosome-Level Assembly

#### The docker image of GABOLA is freely accessible on [**our dockerhub website**](https://hub.docker.com/r/lsbnb/gabola)

## Usage:

The following are the minimum required commands to run each main module separately:

#### 1. LAB Gap Filling:
```
  python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
  /opt/LAB_GapFilling/ProduceBXList.sh -f draft_bwa_mem_C70M60.sam -a draft.fa -o LAB_GapFilling/
  /opt/LAB_GapFilling/Assemble.sh -r /10x_preprocess/nonDupFq/split -o LAB_GapFilling/
  /opt/LAB_GapFilling/Fill.sh -a draft.fa -o LAB_GapFilling/
```
#### 2. GCB Gap Filling:
It can be a stand alone module without any 10x Genomics resource.
```
  /opt/GCB_GapFilling/Fill.sh -a draft.fa -x g-contigs.fa -o GCB_GapFilling/

```
#### 3. LAB Scaffolding:
```
  #rename draft.fa by length first to fit our format
  /opt/Rename_byLength.sh -a draft.fa -o LAB_Scaffolding/
  python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_rename.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
  /opt/LAB_Scaffolding/CandidatePair.sh -f draft_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o LAB_Scaffolding/
  /opt/LAB_Scaffolding/Assemble.sh -f draft_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -r 10x_preprocess/nonDupFq/split -o LAB_Scaffolding/
  /opt/LAB_Scaffolding/Scaffolding.sh -f draft_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_rename.fa -o LAB_Scaffolding/
```

#### 4. GCB Scaffolding:
```
  #rename draft.fa by length first to fit our format
  /opt/Rename_byLength.sh -a draft.fa -o GCB_Scaffolding/
  python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_rename.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/ 
  /opt/GCB_Scaffolding/CandidatePair.sh -f draft_rename_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_rename_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o GCB_Scaffolding/
  /opt/GCB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_rename.fa -x g-contigs.fa -o GCB_Scaffolding/
```

The four main modules can be used in any order and iterated several times. For advanced usage see [**TECHNOTES.md**](TECHNOTES.md)


## Recommended Pipelines & Example Usage:
In terms of Gap Filling, we recommend applying GCB before LAB to fill in bigger gaps first.

As for Scaffolding, GCB is preferably used after LAB as a polishing tool.
### I. 10x Genomics Pipeline
If you only have 10x Genomics linked reads at hand, we propose this pipeline

![alt text](https://eln.iis.sinica.edu.tw/lims/files/users/ccshaney/gabola-gabola_10x_pipeline_0708.jpg)

```
       # GCB Gap-Filling
       ### Under the premise that only 10x reads are available, we would use the initial draft assembly as g-contigs
       /opt/GCB_GapFilling/Fill.sh -a draft.fa -x draft.fa -o GCB_GapFilling/
 
       #LAB Gap-Filling
       python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
       /opt/LAB_GapFilling/ProduceBXList.sh -f draft_gcbgf_bwa_mem_C70M60.sam -a draft_gcbgf.fa -o LAB_GapFilling/
       /opt/LAB_GapFilling/Assemble.sh -r 10x_preprocess/nonDupFq/split/ -o LAB_GapFilling/
       /opt/LAB_GapFilling/Fill.sh -a draft_gcbgf.fa -o LAB_GapFilling/
       
       #LAB Scaffolding
       python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf_labgf.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
       /opt/LAB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o LAB_Scaffolding/
       /opt/LAB_Scaffolding/Assemble.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -r 10x_preprocess/nonDupFq/split -o LAB_Scaffolding/
       /opt/LAB_Scaffolding/Scaffolding.s -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf.fa -o LAB_Scaffolding/

       # GCB Scaffolding
       ###Under the premise that only 10x reads are available, we would use the latest version of the draft assembly as G-contigs. You could also use draft assemblies from any stage of the pipeline.
       python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf_labgf_labs_rename.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
       /opt/GCB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o GCB_Scaffolding/
       /opt/GCB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf_labs_rename.fa -x draft_gcbgf_labgf_labs_rename.fa -o GCB_Scaffolding/

```

### II. 10x Genomics & TGS Pipeline

If 10x Genomics linked reads and Third Generation Sequencing long reads are both obtainable, then we suggest this pipeline:
![alt text](https://eln.iis.sinica.edu.tw/lims/files/users/ccshaney/gabola-gabola_tgs_pipeline_0708.jpg)

```
       #GCB Gap-Filling 
       ### If draft.fa is Supernova draft assembly, we recommend you use TGS long reads or TGS draft assembly as g-contigs.fa
       ### If draft.fa is TGS draft assembly, we recommend you use Supernova draft assembly or unused TGS long reads as g-contigs.fa
       /opt/GCB_GapFilling/Fill.sh -a draft.fa -x g-contigs.fa -o GCB_GapFilling/ 

       #LAB Gap-Filling
       python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
       /opt/LAB_GapFilling/ProduceBXList.sh -f draft_gcbgf_bwa_mem_C70M60.sam -a draft_gcbgf.fa -o LAB_GapFilling/
       /opt/LAB_GapFilling/Assemble.sh -r 10x_preprocess/nonDupFq/split -o LAB_GapFilling/
       /opt/LAB_GapFilling/Fill.sh -a draft_gcbgf.fa -o LAB_GapFilling/

       #LAB Scaffolding
       python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf_labgf.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/ 
       /opt/LAB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o LAB_Scaffolding/
       /opt/LAB_Scaffolding/Assemble.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -r 10x_preprocess/nonDupFq/split -o LAB_Scaffolding/
       /opt/LAB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf.fa -o LAB_Scaffolding/

       #GCB Scaffolding
       ### If draft.fa is Supernova draft assembly, we recommend you use TGS long reads or scaffolds from TGS draft assembly that weren’t used in GCB Gap Filling as g-contigs.fa
       ### If draft.fa is TGS draft assembly, we recommend you use Supernova draft assembly or unused TGS long reads as g-contigs.fa
       python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf_labgf_labs_rename.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/ 
       /opt/GCB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o GCB_Scaffolding/
       /opt/GCB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf_labs_rename.fa -x g-contigs.fa -o GCB_Scaffolding/
```

### III. 10x Genomics & Bionano Genomics Pipeline

Pipeline incorporating 10x Genomics linked reads and Bionano Hybrid Scaffold:
![alt text](https://eln.iis.sinica.edu.tw/lims/files/users/ccshaney/gabola-gabola_bionano_pipeline_0708.jpg)

```
    # GCB Gap-Filling
    ### We recommend you take scaffolds that are in conflict with and those unused by the Bionano protocol (Hybrid Scaffold only utilizes scaffolds longer than 100kbp) to fill in larger gaps in the Hybrid Scaffold first.
    /opt/GCB_GapFilling/Fill.sh -a HybridScaffold.fa -x conflict_unused.fa -o GCB_GapFilling/ 

    # LAB Gap-Filling
    ### Before running this module, scaffolds in conflict and unused by Bionano should be retrieved; therefore, draft_gcbgf.fa should consist of the Hybrid Scaffold after XCB Gap Filling and conflict_unused.fa in XCB Gap Filling
    python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
    /opt/LAB_GapFilling/ProduceBXList.sh -f draft_gcbgf_bwa_mem_C70M60.sam -a draft_gcbgf.fa -o LAB_GapFilling/
    /opt/LAB_GapFilling/Assemble.sh -r 10x_preprocess/nonDupFq/split -o LAB_GapFilling/
    /opt/LAB_GapFilling/Fill.sh -a draft_gcbgf.fa -o LAB_GapFilling/

    #LAB Scaffolding
    python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf_labgf.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/ 
    /opt/LAB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o LAB_Scaffolding/
    /opt/LAB_Scaffolding/Assemble.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -r 10x_preprocess/nonDupFq/split -o LAB_Scaffolding/
    /opt/LAB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf.fa -o LAB_Scaffolding/

    #GCB Scaffolding
    python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_xcbgf_labgf_labs_rename.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/ 
    /opt/XCB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o XCB_Scaffolding/
    ### G-contigs.fa could be the draft assembly or scaffolds that are in conflict with and those unused by the Bionano protocol
    /opt/XCB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf_labs_rename.fa -x xeno-contigs.fa -o GCB_Scaffolding/
```
### IV. 10x Genomics & Bionano Genomics & TGS Pipeline
With all resources mentioned above (10x Genomics linked reads, TGS reads and Bionano Hybrid Scaffold) available, we suggest this pipeline
![alt text](https://eln.iis.sinica.edu.tw/lims/files/users/ccshaney/gabola-gabola_tgs_bionano_pipeline_0708.jpg)

```
     #GCB Gap-Filling
     ### G-contigs.fa could contain (1) Scaffolds in conflict with Bionano cmap (2) Scaffolds unused by Bionano (3) TGS long reads or draft assembly
     We recommend you fill in larger gaps on the Hybrid Scaffold with g-contigs.fa first.
     /opt/GCB_GapFilling/Fill.sh -a HybridScaffold.fa -x g-contigs.fa -o GCB_GapFilling/ 

     #LAB Gap-Filling
     ### Before running this module, scaffolds in conflict and unused by Bionano should be retrieved; therefore, draft_gcbgf.fa should consist of the Hybrid Scaffold after GCB Gap Filling and (1) & (2) aforementioned
     python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
     /opt/LAB_GapFilling/ProduceBXList.sh -f draft_gcbgf_bwa_mem_C70M60.sam -a draft_gcbgf.fa -o LAB_GapFilling/
     /opt/LAB_GapFilling/Assemble.sh -r 10x_preprocess/nonDupFq/split -o LAB_GapFilling/
     /opt/LAB_GapFilling/Fill.sh -a draft_gcbgf.fa -o LAB_GapFilling/

     #LAB Scaffolding
     python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_gcbgf_labgf.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
     /opt/LAB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p  draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o LAB_Scaffolding/
     /opt/LAB_Scaffolding/Assemble.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -r 10x_preprocess/nonDupFq/split -o LAB_Scaffolding/
     /opt/LAB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf.fa -o LAB_Scaffolding/

     #GCB Scaffolding
     python /opt/10x_program/runStep1to3.py -f raw_fastq_dir/ -g draft_xcbgf_labgf_labs_rename.fa --id PROJECTID -a bwa_mem -o 10x_preprocess/
     /opt/GCB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o XCB_Scaffolding/
     ### Since unused and conflict scaffolds were already added back to the draft assembly. G-contigs.fa could be TGS long reads or scaffolds of draft assembly that weren’t used in GCB Gap Filling
     /opt/GCB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf_labs_rename.fa -x g-contigs.fa -o GCB_Scaffolding/
```
