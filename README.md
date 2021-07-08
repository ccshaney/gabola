# GABOLA: A Reliable Gap-Filling Strategy for de novo Chromosome-Level Assembly

##### The docker image of GABOLA is freely accessible at [**our dockerhub website**](https://hub.docker.com/r/lsbnb/gabola)

## Usage:

## Recommended Pipelines & Example Usage:
### I. 10x Genomics Pipeline
If you only have 10x Genomics linked reads at hand, we propose this pipeline

![alt text](https://eln.iis.sinica.edu.tw/lims/files/users/ccshaney/gabola-gabola_10x_pipeline_0707.jpg)

```
       python /opt/10x_program/step1_preprocessFastq.py -f 10x_rawread.fastq --id PROJECTID -o 10x_preprocess/

       #GCB Gap-Filling
       ###Under the premise that only 10x reads are available, we would use the initial draft assembly as xeno-contigs
       /opt/GCB_GapFilling/Fill.sh -a draft.fa -x draft.fa -o GCB_GapFilling/
 
       #LAB Gap-Filling
       python /opt/10x_program/step2_alignment_andFilter.py -a bwa_mem -g draft_gcbgf.fa -f1 NonDupR1.fq -f2 NonDupR2.fq -o draft_gcbgf
       /opt/LAB_GapFilling/ProduceBXList.sh -f draft_gcbgf_bwa_mem_C70M60.sam -a draft_gcbgf.fa -o LAB_GapFilling/
       /opt/LAB_GapFilling/Assemble.sh -r 10x_preprocess/nonDupFq/split -o LAB_GapFilling/
       /opt/LAB_GapFilling/Fill.sh -a draft_gcbgf.fa -o LAB_GapFilling/
       
       #LAB Scaffolding
       python /opt/10x_program/step2_alignment_andFilter.py -a bwa_mem -g draft_gcbgf_labgf.fa -f1 NonDupR1.fq -f2 NonDupR2.fq -o draft_gcbgf_labgf
       python /opt/10x_program/step3_process_samfile.py -f draft_gcbgf_labgf.fa -s draft_gcbgf_labgf_bwa_mem_C70M60.sam 
       /opt/LAB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o LAB_Scaffolding/
       /opt/LAB_Scaffolding/Assemble.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -r 10x_preprocess/nonDupFq/split -o LAB_Scaffolding/
       /opt/LAB_Scaffolding/Scaffolding.s -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf.fa -o LAB_Scaffolding/

       #GCB Scaffolding
       python /opt/10x_program/step2_alignment_andFilter.py -a bwa_mem -g draft_gcbgf_labgf_labs_rename.fa -f1 NonDupR1.fq -f2 NonDupR2.fq -o draft_gcbgf_labgf_labs
       python /opt/10x_program/step3_process_samfile.py -f draft_gcbgf_labgf_labs_rename.fa -s draft_gcbgf_labgf_labs_bwa_mem_C70M60.sam 
       /opt/GCB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o GCB_Scaffolding/
       /opt/GCB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf_labs_rename.fa -x draft_gcbgf_labgf_labs_rename.fa -o GCB_Scaffolding/
       ###Under the premise that only 10x reads are available, we would use the latest version of the draft assembly as G-contigs. You could also use draft assemblies from any stage of the pipeline.

```

### II. 10x Genomics & TGS Pipeline

If 10x Genomics linked reads and Third Generation Sequencing long reads are both obtainable, then we suggest this pipeline:
![alt text](https://eln.iis.sinica.edu.tw/lims/files/users/ccshaney/gabola-gabola_tgs_pipeline_0707.jpg)

```
       python /opt/10x_program/step1_preprocessFastq.py -f 10x_rawread.fastq --id PROJECTID -o 10x_preprocess/

       #GCB Gap-Filling 
       ### If draft.fa is Supernova draft assembly, we recommend you use TGS long reads or TGS draft assembly as g-contigs.fa
       ### If draft.fa is TGS draft assembly, we recommend you use Supernova draft assembly or unused TGS long reads as g-contigs.fa
       /opt/GCB_GapFilling/Fill.sh -a draft.fa -x g-contigs.fa -o XCB_GapFilling/ 

       #LAB Gap-Filling
       python /opt/10x_program/step2_alignment_andFilter.py -a bwa_mem -g draft_gcbgf.fa -f1 NonDupR1.fq -f2 NonDupR2.fq -o draft_gcbgf
       /opt/LAB_GapFilling/ProduceBXList.sh -f draft_gcbgf_bwa_mem_C70M60.sam -a draft_gcbgf.fa -o LAB_GapFilling/
       /opt/LAB_GapFilling/Assemble.sh -r 10x_preprocess/nonDupFq/split -o LAB_GapFilling/
       /opt/LAB_GapFilling/Fill.sh -a draft_gcbgf.fa -o LAB_GapFilling/

       #LAB Scaffolding
       python /opt/10x_program/step2_alignment_andFilter.py -a bwa_mem -g draft_gcbgf_labgf.fa -f1 NonDupR1.fq -f2 NonDupR2.fq -o draft_xcbgf_labgf
       python /opt/10x_program/step3_process_samfile.py -f draft_gcbgf_labgf.fa -s draft_gcbgf_labgf_bwa_mem_C70M60.sam 
       /opt/LAB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o LAB_Scaffolding/
       /opt/LAB_Scaffolding/Assemble.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -r 10x_preprocess/nonDupFq/split -o LAB_Scaffolding/
       /opt/LAB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf.fa -o LAB_Scaffolding/

       #GCB Scaffolding
       python /opt/10x_program/step2_alignment_andFilter.py -a bwa_mem -g draft_gcbgf_labgf_labs_rename.fa -f1 NonDupR1.fq -f2 NonDupR2.fq -o draft_gcbgf_labgf_labs
       python /opt/10x_program/step3_process_samfile.py -f draft_gcbgf_labgf_labs_rename.fa -s draft_gcbgf_labgf_labs_bwa_mem_C70M60.sam 
       /opt/GCB_Scaffolding/CandidatePair.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt.tsv -p draft_gcbgf_labgf_bwa_mem_C70M60_ScafHeadTail_BX_pairSum.tsv -o GCB_Scaffolding/
       ### If draft.fa is Supernova draft assembly, we recommend you use TGS long reads or scaffolds from TGS draft assembly that weren’t used in GCB Gap Filling as g-contigs.fa
       ### If draft.fa is TGS draft assembly, we recommend you use Supernova draft assembly or unused TGS long reads as g-contigs.fa
       /opt/GCB_Scaffolding/Scaffolding.sh -f draft_gcbgf_labgf_labs_bwa_mem_C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv -a draft_gcbgf_labgf_labs_rename.fa -x g-contigs.fa -o GCB_Scaffolding/
```
