#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
This generates a report of variants matching Clinvar variants in 384 genes with zygosity and compound hets annotated
"""
import tiledb.cloud
import math
import pandas as pd
import numpy
import tiledb
from tiledb.cloud.dag import dag
import os
from datetime import datetime
from cryptography.fernet import Fernet
import pyarrow
from typing import Iterable, List, Optional
import logging
logger = logging.getLogger(__name__)

pd.set_option('display.max_rows', None)

#TileDB-VCF array to use
TILEDB_VCF_URI = {}
TILEDB_VCF_URI['RADYPATIENTS'] = "tiledb://Rady_Childrens_Institute_Genomic_Medicine/####-####-####-####-####"
TILEDB_VCF_URI['PRELIMINARY-PHASE-2-PROSPECTIVE-SAMPLES-RE-PROCESSING'] = "tiledb://Rady_Childrens_Institute_Genomic_Medicine/####-####-####-####-####"
VEP_URI = "tiledb://Rady_Childrens_Institute_Genomic_Medicine/vep-variants-dragen-3-9-3-hg38-graph-based-vcfs"
SAMPLE_METADATA_URI = 'tiledb://Rady_Childrens_Institute_Genomic_Medicine/####-####-####-####-####'

GNOMAD_3_1_URI = "tiledb://Rady_Childrens_Institute_Genomic_Medicine/gnomad-annotations"
GNOMAD_4_0_URI = "tiledb://TileDB-Inc/gnomad-4_0-include-nopass"
GNOMAD_URI = GNOMAD_3_1_URI

#use this variant file
VARIANT_SELECTION = 'ANNOTATED_VARIANT_DB'
#VCF selection
VCF_SELECTION = 'RADYPATIENTS'

#the csv variants of interest file in tiledb
VARIANT_FILE_URI = {}

VARIANT_FILE_URI['ANNOTATED_VARIANT_DB'] = "../../data/joint_variants.csv"
VARIANT_FILE = VARIANT_FILE_URI[VARIANT_SELECTION] 

#where the filedb array that holds the variant list is stored
VARIANT_URI_BASE_NAME = {}

# update this if the variant list OR THE BLOCKLIST is altered
VARIANT_URI_BASE_NAME['ANNOTATED_VARIANT_DB'] = "ANNOTATED_VARIANT_DB"         #
VARIANT_URI = 'tiledb://Rady_Childrens_Institute_Genomic_Medicine/s3://tiledb-rchsd-arrays/'+VARIANT_URI_BASE_NAME[VARIANT_SELECTION]

#where the report is stored
adjust_position_for_clinvar = False

#GQ/DP filters
FILTER_FOR_QUALITY = False

FINAL_REPORT_BASE_NAME = {}

FINAL_REPORT_BASE_NAME['ANNOTATED_VARIANT_DB'] = "ANNOTATED_VARIANT_DBfilteredreport"

import datetime
current_datetime = datetime.datetime.now()
formatted_datetime = current_datetime.strftime("%Y%m%d%H%M%S")
REPORT_URI = 'tiledb://Rady_Childrens_Institute_Genomic_Medicine/s3://tiledb-rchsd-arrays/production/'+FINAL_REPORT_BASE_NAME[VARIANT_SELECTION]+"_"+formatted_datetime

#patterns of inheritance
NBS_POI_FILE = "../../data/phase2_beginngs_moi.txt"

BLOCKLIST_FILE = "../../data/prelim_blocklist_dups_removed.csv"

USE_BLOCKLIST = False

SAMPLE_USE_FILE = '../../data/rady_sample_approved_addl_use_indicator.csv'

today_date = datetime.datetime.now().strftime("%d%m%Y")
namespace = "Rady_Childrens_Institute_Genomic_Medicine"

#convert to standard hugo gene names
convert_to_hugo = False
debug = False

#turn off future warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# In[ ]:


#create a dataframe with name, description, stage, unit, value
metric = pd.DataFrame(columns=['name','description','stage','unit','value'])


# In[ ]:


DAG_MODE = tiledb.cloud.dag.Mode.REALTIME

#these are only useful for batch mode
DEFAULT_RESOURCES = {"cpu": "8", "memory": "4Gi"}
DEFAULT_IMG_NAME = "3.9-genomics"


# In[ ]:


def query_sample_metadata(
    sample_stems: list,
    family_roles: list,
    genome_id_list: pyarrow.lib.StringArray,
    remove_sample_name_suffix: bool = True,
) -> pyarrow.Table:
    """
    This function queries the metadata to lookup sample names from genome ids
    Arguments:
        @param genome_id_list: the genome_ids of interest
    Returns:
        Sample data subset
    """
    import pandas
    import pyarrow
    import tiledb

    sample_metadata_attrs = [
        "site_id",
        "study_id",
        "case_id",
        "individual_id",
        "biospecimen",
        "aliquot_id",
        "test_type",
        "family_relationship",
        "version",
        "sample_name",
        "gender",
        "fabric_report_id",
    ]
    with tiledb.open(SAMPLE_METADATA_URI) as A:
        df = A.query(attrs=sample_metadata_attrs)[:]
        df = pandas.DataFrame(df)
        if sample_stems is not None and len(sample_stems) > 0:
            df = df[df.sample_name.str.contains("|".join(sample_stems))]
        if family_roles is not None and len(family_roles) > 0:
            df = df[df.family_relationship.isin(family_roles)]
        if genome_id_list is not None and len(genome_id_list) > 0:
            genome_ids = genome_id_list.to_pylist()
            df = df[df.genome_id.isin(genome_ids)]
        #remove -1 from end of sample names
        if remove_sample_name_suffix:
            df["sample_name"] = df["sample_name"].str.replace("-1$", "")#, regex=True)
    return pyarrow.Table.from_pandas(df)


# This fetches the subset of samples we are interested in. It also fetches sample metadata.
# Because the sample metadata has been comprimised we are not going to use this at the moment as a sample subset.

# In[ ]:


sample_graph = dag.DAG(max_workers=1, namespace=namespace, name = "fetch samples nbs bill")
filtered_metadata_node = sample_graph.submit(
    query_sample_metadata,
    sample_stems=None,
    family_roles=None,
    genome_id_list=None,
    remove_sample_name_suffix=True,
    name=f"Sample Metadata Filter for Probands",
    resource_class="standard",
    
)
sample_graph.compute()
sample_graph.wait()
filtered_metadata_pre_use = filtered_metadata_node.result().to_pandas()


# In[ ]:


#metric = pd.DataFrame(columns=['name','description','stage','unit','value'])
#add the filtered metadata count
metric = metric.append({'name':'filtered_metadata_count','description':'Number of samples in filtered metadata','stage':'pre-use','unit':'count','value':filtered_metadata_pre_use.shape[0]},ignore_index=True)


# Merge with usage indicator file

# In[ ]:


sample_use = pd.read_csv(SAMPLE_USE_FILE)


# In[ ]:


filtered_metadata_use = filtered_metadata_pre_use.merge(sample_use, how='inner', left_on=['individual_id','family_relationship','case_id'], right_on=['ind_id','family_relationship','case_id'])


# In[ ]:


metric = metric.append({'name':'filtered_metadata_count','description':'Number of samples in filtered metadata','stage':'post-sample-use','unit':'count','value':filtered_metadata_use.shape[0]},ignore_index=True)


# In[ ]:


filtered_metadata = filtered_metadata_use[filtered_metadata_use['additional_use'] == 'Yes']


# In[ ]:


metric = metric.append({'name':'filtered_metadata_count','description':'Number of samples in filtered metadata','stage':'additional_use is Yes','unit':'count','value':filtered_metadata.shape[0]},ignore_index=True)

Let's make sure there are no -1's in these sample names 
# In[ ]:


sm_samples = list(filtered_metadata['sample_name'])
len(sm_samples)


# In[ ]:


if not os.path.isfile(NBS_POI_FILE):
    tiledb.cloud.file.export_file("tiledb://Rady_Childrens_Institute_Genomic_Medicine/####-####-####-####",output_uri=NBS_POI_FILE)
else:
    print(f"NBS_POI_FILE already exists {NBS_POI_FILE}")


# This read_partition is meant to be run on a gene basis at the lowest granularity, so it can find compound hets

# In[ ]:


def read_partition(
        uri: str,
        regions_by_gene: List[List[str]],
        sample_partition_count: int,
        sample_partition_id: int,
        *,
        samples: Optional[List[str]] = None,
        region_partition_count: Optional[int] = None,
        selected_variants_uri: Optional[str] = None,
        region_partition_id: Optional[int] = None,
        multiallelic_treatment: str = "decompose",
        drop_duplicates: bool = True,
        remove_sample_name_suffix=True,
        subtract_start_pos_by_one_for_deletions: bool = False,
        variant_selection = None,
        filter_for_quality = False,
        **kwargs
) -> pyarrow.Table:
    """Reads a single partition's worth of VCF data.

    multiallelic_treatment: "hide" or "decompose"
    """
    import pprint
    import tiledbvcf
    import pandas as pd
    import pyarrow as pa
    import re


    cfg = tiledbvcf.ReadConfig(
        sample_partition=(sample_partition_id, sample_partition_count)
    )

            
    # Printed output will be available in the task logs within TileDB cloud
    # once the task is complete. (It will not be seen locally.)
    print("loading regions:")
    pprint.pprint(f"{region_partition_id} of {region_partition_count}")
    print("with configuration:")
    pprint.pprint(cfg)

    vcf_ds = tiledbvcf.Dataset(uri, mode="r", stats=True, cfg=cfg)

    #attrs = vcf_ds.attributes()
    attrs = ["sample_name", "contig", "pos_start", "pos_end", "alleles", "fmt_GT", "fmt_DP","fmt_GQ"]
    if multiallelic_treatment == 'hide':
        attrs += ["info_AN"]
    #this is mostly boilerplate we use for validating input ranges
    clean_regions = []
    contigs = []
    
    #remember regions_by_gene is a list snp by a list of genes
    #if you really have more than one region do the gene, else do all genes
    if region_partition_count is not None and region_partition_count>1:
        
        # break the regions into partitions
        regions_by_partitions = numpy.array_split(regions_by_gene, region_partition_count)
        region_in_partition = regions_by_partitions[region_partition_id]
        regions = [snp for gene in region_in_partition for snp in gene]
    else:
        regions = [snp for gene in regions_by_gene for snp in gene]

    for region in regions:
        regexgroups = re.match("(chr[0-9XYMT]+):([0-9]+)(-([0-9]+))?", region.replace(" ", ""))
        if regexgroups is None:
            print("Invalid genomic coordinate: {genomic_coodinate}")
            return None
        regexchr = regexgroups.group(1)
        regexstart = int(regexgroups.group(2))
        if regexgroups.group(3) is None:
            regexend = regexstart  # point coordinate
        else:
            regexend = int(regexgroups.group(4))  # lose the dash
        contigs += [regexchr]
        clean_regions += [f"{regexchr}:{regexstart}-{regexend}"]
    
    clean_string = ','.join(clean_regions)[0:15]
    print(f"loading {clean_string}...{len(clean_regions)} regions")
    vcf_dfs = []
    vcf_dfs.append(vcf_ds.read_arrow(
        attrs=attrs,
        regions=clean_regions,
        samples=samples
    ))
    while not vcf_ds.read_completed():
        print("still more to go...")
        vcf_dfs.append(vcf_ds.continue_read_arrow())
    
    vcf_df = pyarrow.concat_tables(vcf_dfs).to_pandas()
    

    if vcf_df is None or len(vcf_df) == 0:
        return pyarrow.Table.from_pandas(pandas.DataFrame())
    else:
        print(f"{len(vcf_df)} results in vcf")
    
    alleles = pd.DataFrame(vcf_df["alleles"].to_list())
    vcf_df["ref"] = alleles[0]
    vcf_df["alt"] = alleles.drop([0], axis=1).values.tolist()
    
    if remove_sample_name_suffix == True:
        vcf_df["sample_name"] = vcf_df["sample_name"].str.replace("-1$", "") 
    
    #filter for depth and genotype quality
    if filter_for_quality == True:
        vcf_df = vcf_df[(vcf_df['fmt_GQ'] > 20) & (vcf_df['fmt_DP'] > 8)]
    
    if multiallelic_treatment == 'hide':
        vcf_df = vcf_df.loc[vcf_df['info_AN'] < 3]
    
    # multiple alleles get their own rows
    if multiallelic_treatment == 'decompose':
        vcf_df = vcf_df.explode("alt")
        vcf_df = vcf_df[vcf_df["alt"].notnull()]

    if drop_duplicates:
        #make these "hashable" for downstream operations
        vcf_df["alleles"] = vcf_df["alleles"].apply(lambda x: tuple(x))
        vcf_df["fmt_GT"] = vcf_df["fmt_GT"].apply(lambda x: tuple(x))
        vcf_df = vcf_df.drop_duplicates()
        print(f"{len(vcf_df)} results in vcf after dedup")

    def classifyZygosity(x):
        #rady did not ingest gVCFs so HOM_REF is not possible
        if tuple(x) == tuple([0, 0]):
            raise Exception("HOM_REF not possible")
        elif tuple(x) == tuple([0, 1]):
            return "HET"
        elif tuple(x) == tuple([1, 1]):
            return "HOM_ALT"
        elif tuple(x) == tuple([1]):
            return "HEMI"
        else:
            if multiallelic_treatment == 'hide':
                raise Exception(f"UNKNOWN zygosity {x}")
            else:
                if x[0]==x[1]:
                    return "HOM_ALT"
                else:
                    #look for other edge cases
                    return "HET"
    vcf_df['zygosity'] = vcf_df['fmt_GT'].apply(classifyZygosity)

    #we still need to match REF/ALT
    config = tiledb.Config()
    # 11.53 MB buffer size, is this the culprit?
    config["py.init_buffer_bytes"] = 1024**2 * 11
    ctx = tiledb.Ctx(config=config)

    if selected_variants_uri is not None:
        with tiledb.open(selected_variants_uri,ctx=ctx) as selected_variants:
            sv_df_res = selected_variants.df[:]
            #clinvar pickle has these columns
            if variant_selection == 'CLINVARFILTERED':
                if(set(['contig','pos_start'])).issubset(sv_df_res.columns):
                    sv_df_res = sv_df_res.rename(columns={'contig':'CHROM','pos_start':'POS','GENEINFO':'GENE','pos_end':'clin_end'})
            elif variant_selection == 'ANNOTATED_VARIANT_DB':
                #select these columns: CHROM POS GENE clin_end CLINVAR_VARIANT_ID
                sv_df_res = sv_df_res[['CHROM','POS', 'REF','ALT','GENE','CLINVAR_ID']]
            
            assert(set(['CHROM','POS','REF','ALT']).issubset(sv_df_res.columns))
            assert(set(['contig','pos_start','ref','alt']).issubset(vcf_df.columns))

            if subtract_start_pos_by_one_for_deletions:
                sv_df_res['ref_len'] = sv_df_res['REF'].str.len()
                sv_df_res['alt_len'] = sv_df_res['ALT'].str.len()
                sv_df_res.rename(columns={'POS':'POS_UNADJUSTED'},inplace=True)
                def subtract_start_pos(row):
                    if row['ref_len'] > row['alt_len']:
                        return row['POS_UNADJUSTED'] - 1
                    else:
                        return row['POS_UNADJUSTED']

                # Apply the lambda function to the DataFrame
                sv_df_res['POS'] = sv_df_res.apply(subtract_start_pos, axis=1)


            if len(sv_df_res) == 0:
                return emptyPaTable
            else:
                print(f"{len(sv_df_res)} results in the selected variants")

            sv_df_res['POS']=sv_df_res['POS'].astype(int)
            vcf_df['pos_start']=vcf_df['pos_start'].astype(int)
            vcf_df = vcf_df.merge(sv_df_res,left_on=['contig','pos_start','ref','alt'],right_on=['CHROM','POS','REF','ALT'],how="inner")
            
            if len(vcf_df) == 0:
                print("no selected variants in the vcf")
                return emptyPaTable
            else:
                print(f"{len(vcf_df)} variant-sample tuples in the vcf after merge with selected variants")
            
            hets = vcf_df[vcf_df['zygosity'] == 'HET']

            
            if len(hets) > 0:
                grouped_hets = hets.groupby(['sample_name','GENE']).size().reset_index(name='het_count')
                compound_hets = grouped_hets[grouped_hets['het_count'] > 1].copy()
                if len(compound_hets) > 0:
                    #avoid the A value is trying to be set on a copy of a slice from a DataFrame error
                    compound_hets.loc[:,'compound_event']='CMPD_HET'
                    compound_hets.loc[:,'zygosity']='HET' #for the join
                    print(f"{len(compound_hets)} compound hets found")
                    vcf_df = vcf_df.merge(compound_hets,on=['sample_name','zygosity','GENE'],how="left")
                else:
                    #columns should match up
                    vcf_df.loc[:,'het_count'] = 0
                    vcf_df.loc[:,'compound_event']=''
            else:
                vcf_df.loc[:,'het_count'] = 0
                vcf_df.loc[:,'compound_event']=''
            
            vcf_df['het_count'] = vcf_df['het_count'].fillna(0)
            vcf_df['het_count'] = vcf_df['het_count'].astype('int64')
            
    
    if len(vcf_df) == 0:
        return emptyPaTable
    else:
        print(f"{len(vcf_df)} variant-sample tuples in vcf")
    
    arrow_table = pa.Table.from_pandas(vcf_df)
    assert(isinstance(arrow_table,pyarrow.lib.Table))
    return arrow_table


# In[ ]:


def combine_results(df_list, *args, **kwargs):
    import pyarrow

    print(f"Input list contains {len(df_list)} items")
    results = []
    for x in df_list:
        if x is None:
            continue
        if not isinstance(x,pyarrow.lib.Table):
            #this can happen when pyarrow gets emotionally overwhelmed
            print("Not returing a pyarrow table, instead its a {type(x)}")
            continue
        if x.num_rows>0:
            results+=[x]
        else:
            continue
    print('There are ', len(results), 'non-empty results')
    if len(results) > 1:
        table = pyarrow.concat_tables(results)
    elif len(results) == 1:
        table = results[0]
    else:
        table = None
    return table


# In[ ]:


import pandas as pd
import tiledb
import numpy as np


dtypes_clinvar = {'internal_id':np.int32, 'contig':str, 'pos_start':np.int32, 
       'CLINVAR_ID':str, 'REF':str, 'ALT':str, 'ALLELEID':str,
       'CLNDISDB':str, 'CLNDN':str, 'CLNHGVS':str, 'CLNREVSTAT':str, 'CLNSIG':str, 'CLNVC':str,
       'CLNVCSO':str, 'GENEINFO':str, 'ORIGIN':str, 'MC':str, 'pos':str, 'pos_end':str, 'pos_37':str,'pos_start_37':object,'pos_end_37':object,
       'NC':str, 'G_DOT':str,
       'NM':str, 'C_DOT':str, 'NP':str, 'P_DOT':str,'sub_classification':str}

#Rady.list,inheritance.list,ID,CHROM,POS,REF,ALT,CLINVAR_VARIANT_ID,CLINVAR_ALLELE_ID,CLINVAR_DISEASE,CLINVAR_REVIEW_STATUS,CLINVAR_INTERP,GENE,CLINVAR_CONFLICT_CLASS,CLINVAR_INTERP_SIMPLE,MM_INTERP,MM_DISEASE,MM_VARIANT_ID,MM_INTERP_SIMPLE,JOINT_INTERP_SIMPLE,internal_id,CLNHGVS,CLNSIG,CLNVC,pos_38,pos_37,NC,G_DOT,NM,C_DOT,NP,P_DOT,sub_classification,Fabric_chr,Fabric_b37_pos,ACE.gene,ACE.consequence,ACE.classification,ACMG.rules.matched,Variant..UKB_EURpe.MCPS,Ref.Allele.MCPS,Alt.Allele.MCPS,CADD_PHRED.MCPS,Function.MCPS,Gene.Name.MCPS,Major.Hom.Ctrl.MCPS,Het.Ctrl.MCPS,Minor.Hom.Ctrl.MCPS,Minor.Hom.Ctrl.Freq.MCPS,Het.Ctrl.Freq.MCPS,Ctrl.Maf.MCPS,Ctrl.HWE_P.MCPS,Gene.ID,Transcript.Stable.Id,Has.CCDS.Transcript,Transcript.codon.change,Transcript.AA.change,.MCPSon.Rank,Gene.Transcript,LOF.Number.Of.Transcripts.In.Gene.Affected,LOF.Percentage.Of.Transcripts.Affected,OMIM_Number,OMIM_Pheno.UKB_EURpe,ClinVar.Disease,ClinVar.Clinical.Significance,ClinVar.PMID,MTR,MTR.FDR,MTR.Centile,REVEL,MGI_Pheno.UKB_EURpes,MGI_Essential,med_QUAL,med_PercAltRead,med_READS_REF,med_READS_ALT,med_GQ,med_FS,med_MQ,med_QD,med_MQRankSum,med_ReadPosRankSum,med_coverage,med_PL,med_DPF,med_R2_5P_bias,med_AF,med_SOR,med_FractionInformativeReads,Variant..UKB_EURpe.UKB_EUR,CADD_PHRED.UKB_EUR,Is.Minor.Ref,Function.UKB_EUR,Major.Hom.Ctrl.UKB_EUR,Het.Ctrl.UKB_EUR,Minor.Hom.Ctrl.UKB_EUR,Minor.Hom.Ctrl.Freq.UKB_EUR,Het.Ctrl.Freq.UKB_EUR,Missing.Ctrl,QC.Fail.Ctrl,Ctrl.Maf.UKB_EUR,Case.HWE_P,Ctrl.HWE_P.UKB_EUR
fields_genomenon = ['Rady.list', 'inheritance.list', 'ID', 'CHROM', 'POS', 'REF', 'ALT', 'CLINVAR_VARIANT_ID', 'CLINVAR_ALLELE_ID', 'CLINVAR_DISEASE', 'CLINVAR_REVIEW_STATUS', 'CLINVAR_INTERP', 'GENE', 'CLINVAR_CONFLICT_CLASS', 'CLINVAR_INTERP_SIMPLE', 'MM_INTERP', 'MM_DISEASE', 'MM_VARIANT_ID', 'MM_INTERP_SIMPLE', 'JOINT_INTERP_SIMPLE', 'internal_id', 'CLNHGVS', 'CLNSIG', 'CLNVC', 'pos_38', 'pos_37', 'NC', 'G_DOT', 'NM', 'C_DOT', 'NP', 'P_DOT', 'sub_classification', 'Fabric_chr', 'Fabric_b37_pos', 'ACE.gene', 'ACE.consequence', 'ACE.classification', 'ACMG.rules.matched', 'Variant..UKB_EURpe.MCPS', 'Ref.Allele.MCPS', 'Alt.Allele.MCPS', 'CADD_PHRED.MCPS', 'Function.MCPS', 'Gene.Name.MCPS', 'Major.Hom.Ctrl.MCPS', 'Het.Ctrl.MCPS', 'Minor.Hom.Ctrl.MCPS', 'Minor.Hom.Ctrl.Freq.MCPS', 'Het.Ctrl.Freq.MCPS', 'Ctrl.Maf.MCPS', 'Ctrl.HWE_P.MCPS', 'Gene.ID', 'Transcript.Stable.Id', 'Has.CCDS.Transcript', 'Transcript.codon.change', 'Transcript.AA.change', '.MCPSon.Rank', 'Gene.Transcript', 'LOF.Number.Of.Transcripts.In.Gene.Affected', 'LOF.Percentage.Of.Transcripts.Affected', 'OMIM_Number', 'OMIM_Pheno.UKB_EURpe', 'ClinVar.Disease', 'ClinVar.Clinical.Significance', 'ClinVar.PMID', 'MTR', 'MTR.FDR', 'MTR.Centile', 'REVEL', 'MGI_Pheno.UKB_EURpes', 'MGI_Essential', 'med_QUAL', 'med_PercAltRead', 'med_READS_REF', 'med_READS_ALT', 'med_GQ', 'med_FS', 'med_MQ', 'med_QD', 'med_MQRankSum', 'med_ReadPosRankSum', 'med_coverage', 'med_PL', 'med_DPF', 'med_R2_5P_bias', 'med_AF', 'med_SOR', 'med_FractionInformativeReads', 'Variant..UKB_EURpe.UKB_EUR', 'CADD_PHRED.UKB_EUR', 'Is.Minor.Ref', 'Function.UKB_EUR', 'Major.Hom.Ctrl.UKB_EUR', 'Het.Ctrl.UKB_EUR', 'Minor.Hom.Ctrl.UKB_EUR', 'Minor.Hom.Ctrl.Freq.UKB_EUR', 'Het.Ctrl.Freq.UKB_EUR', 'Missing.Ctrl', 'QC.Fail.Ctrl', 'Ctrl.Maf.UKB_EUR', 'Case.HWE_P', 'Ctrl.HWE_P.UKB_EUR']
dtypes_genomenon = {field: str for field in fields_genomenon}
dtypes_genomenon['POS'] = int



dtypes = dtypes_genomenon

if VARIANT_FILE_URI[VARIANT_SELECTION].endswith('.gz'):
       variant_df_full = pd.read_csv(VARIANT_FILE_URI[VARIANT_SELECTION], encoding='unicode_escape',dtype=dtypes,compression='gzip',low_memory=False)
else:
       variant_df_full = pd.read_csv(VARIANT_FILE_URI[VARIANT_SELECTION], encoding='unicode_escape',dtype=dtypes,low_memory=False)

if VARIANT_SELECTION == 'GENOMENON' or VARIANT_SELECTION == 'ANNOTATED_VARIANT_DB':
       #rename CLINVAR_VARIANT_ID to CLINVAR_ID
       variant_df_full = variant_df_full.rename(columns={'CLINVAR_VARIANT_ID':'CLINVAR_ID'})

variant_df = variant_df_full[['CHROM','POS','REF','ALT','GENE','CLINVAR_ID']]


# In[ ]:


variant_df_full.head()


# In[ ]:


metric = metric.append({'name':'variants_of_interest','description':'variants in selection list','stage':'initial','unit':'chr-pos-ref-alt','value':variant_df.shape[0]},ignore_index=True)


# In[ ]:


#assert all CHROM starts with chr
assert(variant_df['CHROM'].str.startswith('chr').all())

#assert all pos_start is int
variant_df.loc[:, 'POS'] = variant_df['POS'].astype(int)
assert(variant_df['POS'].dtype == int)

#sort by CHR and POS
variant_df.sort_values(by=['CHROM','POS'],inplace=True)


# In[ ]:


variant_df.head()


# In[ ]:


variant_df=variant_df.replace(np.nan,'',regex=True) 
#show all the types in this dataframe

#which are mixed dtypes?
for col in variant_df.columns:
    weird = (variant_df[[col]].applymap(type) != variant_df[[col]].iloc[0].apply(type)).any(axis=1)
    if len (variant_df[weird]) > 0:
        #set the column to the correct type
        variant_df[col] = variant_df[col].astype(str)


# In[ ]:


variant_df.columns


# In[ ]:


len(variant_df)


# In[ ]:


blocklist=pd.read_csv(BLOCKLIST_FILE)
#rename blocklist columns to CHROM, POS, REF, ALT
blocklist = blocklist.rename(columns={'pos_start':'POS','ref':'REF','alt':'ALT'})
#drop the columns we don't need
blocklist = blocklist[['CHROM','POS','REF','ALT']]
if USE_BLOCKLIST:
       #use CHROM,POS,REF,ALT
       blocked_entries = variant_df.merge(blocklist,how='inner',on=['CHROM','POS','REF','ALT'])

       block_all = variant_df.merge(blocked_entries.drop_duplicates(), on=['CHROM','POS','REF','ALT','GENE','CLINVAR_ID'], 
                   how='left', indicator=True)      
       #remove the entries that are in the blocked_entries
       variant_df = block_all[block_all['_merge'] == 'left_only']


# In[ ]:


len(variant_df)


# In[ ]:


variant_df.columns


# In[ ]:


metric = metric.append({'name':'variants_of_interest','description':'variants in selection list','stage':'after blocklist USE_BLOCKLIST:{}'.format(USE_BLOCKLIST),'unit':'chr-pos-ref-alt','value':variant_df.shape[0]},ignore_index=True)


# In[ ]:



if tiledb.array_exists(VARIANT_URI):
       print("Array already exists, skipping ingest")
else:
       tiledb.from_pandas(uri=VARIANT_URI,dataframe=variant_df,mode="ingest",row_start_idx=0,sparse=False)


# These are to be blocked or masked out of the variant dataset

# In[ ]:


len(variant_df)


# Beware there are some overlapping genes and variants assigned willy-nilly
variant_df[variant_df['POS']==17388025]
# In[ ]:


def regions_snps_of_interest_ranges_by_gene_name(
    gene_name_list: List[str], selected_variants_uri: str
) -> List[List[str]]:
    if gene_name_list is None:
        with tiledb.open(selected_variants_uri) as A:
            sv_ds = A.query(return_arrow=False).df[:]
    else:
        with tiledb.open(selected_variants_uri) as A:
            sv_ds = A.query(return_arrow=False, cond=f"GENE in {gene_name_list}").df[:]
    #get rid of weird artefacts
    sv_df_res_filt = sv_ds[sv_ds['CHROM'].str.startswith('chr')]
    if gene_name_list:
        assert(sorted(list(set(sv_df_res_filt['GENE']))) == sorted(gene_name_list))
    gene_groups = sv_df_res_filt.sort_values(by=["CHROM","POS"]).groupby('GENE')
    regions = [[f'{row["CHROM"]}:{row["POS"]}-{row["POS"]}' for _, row in gene.iterrows()] for _, gene in gene_groups]
    #this table has some entries with different allele changes same position
    for regionidx in range(len(regions)):
        regions[regionidx] = list(set(regions[regionidx]))
    return regions

def regions_snps_of_interest_ranges_by_refseq_name(
    gene_name_list: List[str], selected_variants_uri: str
) -> List[List[str]]:
    if gene_name_list is None:
        with tiledb.open(selected_variants_uri) as A:
            sv_ds = A.query(return_arrow=False).df[:]
    else:
        with tiledb.open(selected_variants_uri) as A:
            sv_ds = A.query(return_arrow=False, cond=f"NM in {gene_name_list}").df[:]
    #get rid of weird artefacts caused by null fill values
    sv_df_res_filt = sv_ds[sv_ds['contig'].str.startswith('chr')]
    #sv_df_res_filt = sv_df_res_filt[~sv_df_res_filt['contig'].str.startswith('chrX')]
    if gene_name_list:
        assert(sorted(list(set(sv_df_res_filt['NM']))) == sorted(gene_name_list))
    gene_groups = sv_df_res_filt.sort_values(by=["contig","pos_start"]).groupby('NM')
    regions = [[f'{row["contig"]}:{row["pos_start"]}-{row["pos_start"]}' for _, row in gene.iterrows()] for _, gene in gene_groups]
    #this table has some entries with different allele changes same position
    for regionidx in range(len(regions)):
        regions[regionidx] = list(set(regions[regionidx]))
    return regions

def wholerange(
        snps_of_interest_ranges: List[List[str]],
):
    allranges = []
    for snps_of_interest_range in snps_of_interest_ranges:
        chrom = snps_of_interest_range[0].split(':')[0]
        # Extract start and end positions from each range
        positions = [list(map(int, r.split(':')[1].split('-'))) for r in snps_of_interest_range]

        # Find the minimum start position and maximum end position
        min_start = min(pos[0] for pos in positions)
        max_end = max(pos[1] for pos in positions)

        # Construct the overall range string
        allranges += [[f'{chrom}:{min_start}-{max_end}']]

    print(allranges)
    return allranges


# Here you can subset the gene list for testing

# In[ ]:


gene_name_list_all = sorted(list(set(variant_df['GENE'])))

gene_name_list = gene_name_list_all

#return a nested list of 1bp regions convering snps per gene from the list of SNPs of interest
if VARIANT_SELECTION == 'CLINVARXML' or VARIANT_SELECTION == 'CLINVARFILTERED':
    regions = regions_snps_of_interest_ranges_by_refseq_name(gene_name_list,VARIANT_URI)
else:
    regions = regions_snps_of_interest_ranges_by_gene_name(gene_name_list,VARIANT_URI)

wholeranges = wholerange(regions)


# In[ ]:


metric = metric.append({'name':'regions_of_interest','description':'regions in selection list by gene or refseq name','stage':'initial','unit':'region','value':len(regions)},ignore_index=True)
metric = metric.append({'name':'whole_ranges_of_interest','description':'whole ranges in selection list min-max','stage':'initial','unit':'whole_range','value':len(wholeranges)},ignore_index=True)


# In[ ]:


#make sure this remains sorted by chromosome
regions = sorted(regions,key=lambda x: x[0].split(':')[0])


# In[ ]:


#save regions to gzipped file, comma separated, one list per line
import gzip
with gzip.open("regions.txt.gz", 'wt') as f:
    for region in regions:
        f.write(','.join(region)+'\n')


# In[ ]:


VARIANT_URI


# In[ ]:


print(f"{len(regions)} genes")


# In[ ]:


print(f"{sum(len(gene) for gene in regions)} loci ")


# Let's use realtime task graph mode for this query

# In[ ]:


graph = tiledb.cloud.dag.DAG(
    max_workers=60,
    namespace=namespace,
    name = "radyANNOTATED_VARIANT_DB",
    mode=DAG_MODE
)


# Fetch and annotate variants on a region (gene) basis, with presumably all samples in each partition
# 
# We just pass regionlist list elements each one all the desired snps in a gene, no need to partition those further so set region_partition_count to None
# 
# Might take 6-10min

# In[ ]:


sample_partition_count = 60
region_partition_count = 1 #1 means pull all the loci at once
all_reads = []
for region_part_id in [0]:
    for sample_part_id in range(sample_partition_count):
        results = graph.submit(
            read_partition,
            uri = TILEDB_VCF_URI[VCF_SELECTION], 
            regions_by_gene = regions,
            samples = None,
            selected_variants_uri = VARIANT_URI,
            region_partition_count = region_partition_count,
            region_partition_id=region_part_id,
            sample_partition_count=sample_partition_count,
            sample_partition_id=sample_part_id,
            drop_duplicates = True,
            subtract_start_pos_by_one_for_deletions = adjust_position_for_clinvar,
            result_format="arrow",
            local_mode=False,
            name=f"Rady samples realtime query (Region {region_part_id}, Sample {sample_part_id})",
            resource_class = "standard",
            mode=DAG_MODE,
            variant_selection=VARIANT_SELECTION,
            filter_for_quality = FILTER_FOR_QUALITY,
        )
        all_reads.append(results)


# In[ ]:


metric = metric.append({'name':'sample_partition_count','description':'technical metadata about task graph partition','stage':'initial','unit':'int','value':sample_partition_count},ignore_index=True)
metric = metric.append({'name':'region_partition_count','description':'technical metadata about task graph partition','stage':'initial','unit':'int','value':region_partition_count},ignore_index=True)


# In[ ]:


metric = metric.append({'name':'regions_of_interest','description':'regions queried in dag','stage':'pre-query','unit':'int','value':len(regions)},ignore_index=True)


# In[ ]:


len(regions)


# Let's not send a list of prespecified variants but just pull down all variants in the region of interest, get VEP'ed and then join using gdot

# In[ ]:


regions


# In[ ]:


for region in regions:
    print(f"{len(region)} snps")


# Combine the results

# In[ ]:


combined_vcf_results = graph.submit(
    combine_results, 
    all_reads, 
    name="Combinesamples", 
    mode=DAG_MODE,
    local_mode=True,
)


# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
#graph.visualize()


# In[ ]:


graph.compute()

gnomad_graph.compute()
# This can take about 10 minutes

# In[ ]:


dev_mode = False
# ... and wait for it to finish.
retry = True
try:
    graph.wait()
except tiledb.cloud.tiledb_cloud_error.TileDBCloudError as tdberror:
    # Allow retry, if not in dev mode
    if not dev_mode:
        logger.error(f"retriable graph error\n{tdberror}")
        retry = True
    else:
        logger.error(f"fatal graph error\n{tdberror}")
        logger.debug(dag.server_logs(graph))

# Retry once
if retry:
    try:
        logger.info("retry graph")
        graph.retry_all()
        graph.wait()
    except tiledb.cloud.tiledb_cloud_error.TileDBCloudError as tdberror:
        logger.error(f"graph error during retry, exiting\n{tdberror}")
        logger.debug(dag.server_logs(graph))


# In[ ]:


#one more retry just in case
graph.retry_all()
graph.wait()


# In[ ]:


if combined_vcf_results.result() is not None:
    vcf_df = combined_vcf_results.result().to_pandas()
else:
    print("No results")


# In[ ]:


vcf_df.head()


# In[ ]:


len(vcf_df)


# In[ ]:


metric = metric.append({'name':'vcf_df','description':'result of the main query','stage':'post-query','unit':'chr-pos-ref-alt-sample no dedup','value':vcf_df.shape[0]},ignore_index=True)

These sample names should also be -1 free
# In[ ]:


vcf_df['sample_name'].unique()


# In[ ]:


vcf_df.columns


# In[ ]:


unique_combinations_df = vcf_df[['contig', 'pos_start']].drop_duplicates()
distinct_count = unique_combinations_df.shape[0]


# In[ ]:


distinct_count


# In[ ]:


len(vcf_df[['sample_name','contig','pos_start','ref','alt']].drop_duplicates())


# In[ ]:


#assert none of the sample-GENE singletons are called cmpd_het
assert(len(vcf_df.groupby(['sample_name','GENE','compound_event']).size().reset_index(name='count').query('count == 1')) == 0)


# In[ ]:


vcf_df.groupby(['sample_name','GENE','compound_event']).size().reset_index(name='count')


# In[ ]:


vcf_df['sample_name'].nunique()


# This should capture 29 ABCC8 sample-variant tuples regardless of other genes in the query

# In[ ]:


if VARIANT_SELECTION=='KINGSMORE':
    ABCC8 = vcf_df[vcf_df['GENE']=='ABCC8']


# In[ ]:


vcf_df.columns


# This is important. If we want to filter for probands, make sure the VCF columns match up with the metadata and set `how = inner`

# In[ ]:


vcf_df = pd.merge(
    vcf_df, filtered_metadata, on=["sample_name"], how="inner"
)


# In[ ]:


filtered_metadata.columns


# In[ ]:


vcf_df['sample_name'].nunique()


# In[ ]:


len(vcf_df[['sample_name','contig','pos_start','ref','alt']].drop_duplicates())


# Actually drop any duplicates

# In[ ]:


vcf_df = vcf_df.drop_duplicates(subset=['sample_name','contig','pos_start','ref','alt'])


# In[ ]:


metric = metric.append({'name':'vcf_df','description':'result of the main query','stage':'post-de-depup','unit':'distinct chr-pos-ref-alt-sample','value':vcf_df.shape[0]},ignore_index=True)


# In[ ]:


import numpy as np
def convert_to_tuple(value):
    if isinstance(value, list):
        return tuple(value)
    elif isinstance(value, np.ndarray):
        return str(value)
    return value

# Iterate over each column and convert lists to tuples
vcf_df_nt = vcf_df.applymap(convert_to_tuple)
vcf_df_nt_nd = vcf_df_nt.drop_duplicates(subset=['sample_name','contig','pos_start','ref','alt'])
len(vcf_df_nt_nd)
# Output the modified DataFrame


# In[ ]:


vcf_df_nt_nd.head()


# In[ ]:


len(vcf_df[['sample_name','contig','pos_start','ref','alt']].drop_duplicates())


# VCFs do not know males have one X chromosome. So set that manually.

# In[ ]:


#indeed chrX looks like 0/1 for males
def set_hemizygous(row):
    if row['gender'] == 'Male' and row['contig'] == 'chrX':
        return 'HEMI'
    else:
        return row['zygosity']

# Apply the custom function to the DataFrame
vcf_df['zygosity'] = vcf_df.apply(set_hemizygous, axis=1)


# In[ ]:


vcf_df['GENE'] = vcf_df['GENE'].apply(lambda x: x.split(':')[0])


# In[ ]:


len(vcf_df)


# In[ ]:


def read_moi_genes():
    nbs_moi_df = pd.read_csv(NBS_POI_FILE, sep='\t')
    nbs_moi_df = nbs_moi_df.rename(columns = {'MOI':'MOI','Gene':'GENE'})
    new_df = nbs_moi_df.drop_duplicates(subset=['GENE'])
    gene_poi_map = {}
    for gene, poi in zip(new_df['GENE'], new_df['MOI']):
        gene_poi_map[gene.upper()] = poi
    return gene_poi_map

gene_poi_map = read_moi_genes()


# In[ ]:


NBS_POI_FILE


# In[ ]:


gene_poi_map


# In[ ]:


vcf_df['moi'] = vcf_df['GENE'].map(gene_poi_map)


# In[ ]:


vcf_df['moi']


# In[ ]:


def positive_genotype(row):
    if row['moi'] == 'Pattern unknown':
        return "Unknown"
    elif row['moi'] == 'AD':
        return "Yes"
    elif row['moi'] == 'AR':
        if row['zygosity'] == 'HOM_ALT' or row['compound_event'] == 'CMPD_HET':
            return "Yes"
        else:
            return "No"
    elif row["moi"] == 'XR':
        if row['gender'] == "Female":
            if row['zygosity'] == 'HOM_ALT' or row['compound_event'] == 'CMPD_HET':
                return "Yes"
        else:
            if row['zygosity'] == 'HEMI':
                return 'Yes'
        return "No"
    elif row["moi"] == 'XD':
        if row['gender'] == "Male":
            if row['zygosity'] == 'HEMI':
                return 'Yes'
        else:
            return "Yes"


vcf_df['positive_genotype'] = vcf_df.apply(positive_genotype, axis=1)


# In[ ]:


metric = metric.append({'name':'vcf_df','description':'vcf_df entries postive wrt moi','stage':'moi_positive','unit':'distinct chr-pos-ref-alt-sample','value':len(vcf_df[vcf_df['positive_genotype']=='Yes'])},ignore_index=True)


# In[ ]:


columns = [col for col in vcf_df.columns if col not in ['individual_id', 'sample_name']]


# In[ ]:


validation_data = pd.DataFrame({
    "genome_id": [949581, 963614, 949581, 963614, 1033213, 1053847, 1033213, 1053847],
    "contig": ["chr18", "chr18","chr18","chr18","chr18","chr18","chr18","chr18"],
    "pos_start": [23568833, 23573531, 23539394, 23544402, 23536736, 23533472, 23539859, 23539402],
    "ref": ["ACT", "GCA", "G", "G", "A", "A", "T", "G"],
    "alt": ["A", "G", "A", "T", "G", "C", "C", "A"],
    "gene_symbol": ["NPC1", "NPC1", "NPC1", "NPC1", "NPC1", "NPC1", "NPC1", "NPC1"],
    "gene_id": ["ENSG00000141458"]*8,
    "canonical_hgvs_c": ["c.451_452del", "c.99_100del", "c.2872C>T", "c.2072C>A", "c.3182T>C", "c.3637T>G", "c.2747A>G", "c.2864C>T"],
    "canonical_hgvs_p": ["p.Ser151PhefsTer18", "p.Ala34IlefsTer23", "p.Arg958Ter", "p.Pro691Gln", "p.Ile1061Thr", "p.Leu1213Val", "p.Asn916Ser", "p.Ser955Phe"],
    "transcrip_id": ["ENST00000269228"]*8,
    "variant_classification": ["P", "P", "P", "LP", "P", "VUS", "VUS", "VUS"]
})


# In[ ]:


import pandas as pd
pd.set_option("display.max_columns", 100)
pd.set_option("display.max_colwidth",500)


# In[ ]:



validation_data.merge(vcf_df,on=['genome_id','contig','pos_start','ref','alt'],how='inner')


# Check against this gold standard set of 300 variants

# In[ ]:


# load the first sheet of the Excel file into a data frame
goldset = pd.read_excel("../../data/restrospective_variants_truth_set_052423.xlsx", sheet_name=0).rename(columns={"start_hg38":"pos_start"})
goldset['contig'] = 'chr' + goldset['chr'].astype(str)
goldset.drop(columns=['chr','sample_name','zygosity','clinvar_variant_id'], inplace=True)
goldset = goldset.astype({'pos_start':'int'})
goldset.drop_duplicates(subset=['contig','pos_start','ref','alt'], inplace=True)
goldset.sort_values(by=['pos_start']).head()


# In[ ]:


goldset.columns


# In[ ]:


chromosome_order = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']


# In[ ]:


goldset['contig'] = pd.Categorical(goldset['contig'], categories=chromosome_order, ordered=True)
goldset['chromosome_order'] = goldset['contig'].map(lambda c: chromosome_order.index(c))

goldmerge = goldset.merge(vcf_df, on=['contig', 'pos_start','ref','alt'], how='left').sort_values(['chromosome_order', 'pos_start'])[['contig','pos_start','ref','alt','gene_symbol','zygosity','sample_name','CLINVAR_ID']]
goldmerge.to_csv('../../results/goldmerge.csv', index=False)
goldmerge.head()


# In[ ]:


goldcount = goldmerge.groupby(['contig','pos_start','ref','alt']).agg(['nunique']).reset_index(drop=False).reset_index()
goldcount.head()


# In[ ]:


metric = metric.append({'name':'vcf_df','description':'goldset entries','stage':'goldset','unit':'distinct chr-pos-ref-alt-sample','value':sum(goldcount['sample_name']['nunique'])},ignore_index=True)


# In[ ]:


hotsamples = int((goldcount['sample_name']>1).sum())
print("{} of these {} variants has at least one sample".format(hotsamples,len(goldcount)))


# Cast the complex types as strings to work with `from_pandas`

# In[ ]:


obj_cols = vcf_df.select_dtypes(include=['object']).columns
vcf_df[obj_cols] = vcf_df[obj_cols].astype(str)


# In[ ]:


len(vcf_df)


# In[ ]:


vcf_df.columns


# In[ ]:


#merge using CHROM,POS,REF,ALT,GENE,CLINVAR_VARIANT_ID,Rady.list
vcf_df = vcf_df.merge(variant_df_full,on=['CHROM','POS','REF','ALT','GENE'],how='inner')


# In[ ]:


metric = metric.append({'name':'vcf_df','description':'result of the main query','stage':'post-anno','unit':'distinct chr-pos-ref-alt-sample','value':vcf_df.shape[0]},ignore_index=True)


# In[ ]:


metric.to_csv('../../results/metrics.tsv',index=False,sep="\t")


# In[ ]:


vcf_df.columns


# In[ ]:


#gzip this file
vcf_df.to_csv('../../results/bulk_moi.tsv.gz', index=False, sep='\t', compression='gzip')


# In[ ]:


vcf_df.to_csv('../../results/bulk_moi.tsv.zip', index=False, sep='\t', compression='zip')


# In[ ]:


# Python program to find MD5 hash value of a file
import hashlib
 
def md5sum(filename):
    md5_hash = hashlib.md5()
    with open(filename,"rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096),b""):
            md5_hash.update(byte_block)
        return(md5_hash.hexdigest())

resultcsvs = []
for f in range(5):
    print(f)
    resultcsvs += [pd.read_csv(f"../../results/{f}_bulk_moi.tsv.gz",sep="\t",compression="gzip",low_memory=False)[['sample_name','contig','pos_start','ref','alt']].drop_duplicates()]
# In[ ]:


#save this info to a metadata file
run_metadata = {
    "VARIANT_SELECTION":VARIANT_SELECTION,
    "VARIANT_FILE_URI":VARIANT_FILE_URI[VARIANT_SELECTION],
    "VCF_SELECTION":VCF_SELECTION,
    "TILEDB_VCF_URI":TILEDB_VCF_URI[VCF_SELECTION],
    "NBS_POI_FILE":NBS_POI_FILE,
    "BLOCKLIST_FILE":BLOCKLIST_FILE,
    "USE_BLOCKLIST":USE_BLOCKLIST,
    "SAMPLE_USE_FILE":SAMPLE_USE_FILE,
    "START_POS_ADJUSTMENT":adjust_position_for_clinvar,
    "FILTER_FOR_QUALITY":FILTER_FOR_QUALITY,
    "bulk_moi":'../../results/bulk_moi.tsv.gz',
    "goldmerge":'../../results/goldmerge.csv'
}
#fetch the md5sum of each file in the run_metadata
run_metadata['VARIANT_FILE_MD5'] = md5sum(VARIANT_FILE_URI[VARIANT_SELECTION])
run_metadata['NBS_POI_FILE_MD5'] = md5sum(NBS_POI_FILE)
run_metadata['BLOCKLIST_FILE_MD5'] = md5sum(BLOCKLIST_FILE)
run_metadata['SAMPLE_USE_FILE_MD5'] = md5sum(SAMPLE_USE_FILE)
run_metadata['bulk_moi_MD5'] = md5sum('../../results/bulk_moi.tsv.gz')
run_metadata['goldmerge_MD5'] = md5sum('../../results/goldmerge.csv')

#get the line counts of each file in the run_metadata
run_metadata['VARIANT_FILE_LINECOUNT'] = !wc -l {VARIANT_FILE_URI[VARIANT_SELECTION]} | cut -d' ' -f1
run_metadata['NBS_POI_FILE_LINECOUNT'] = !wc -l {NBS_POI_FILE} | cut -d' ' -f1
run_metadata['BLOCKLIST_FILE_LINECOUNT'] = !wc -l {BLOCKLIST_FILE} | cut -d' ' -f1
run_metadata['SAMPLE_USE_FILE_LINECOUNT'] = !wc -l {SAMPLE_USE_FILE} | cut -d' ' -f1
run_metadata['bulk_moi_LINECOUNT'] = !gunzip -c ../../results/bulk_moi.tsv.gz | wc -l | cut -d' ' -f1
run_metadata['goldmerge_LINECOUNT'] = !wc -l ../../results/goldmerge.csv | cut -d' ' -f1

import json
json.dump(run_metadata, open('../../results/run_metadata.json', 'w'), indent=4)

import tiledb
#turn BLOCK_LIST_2022 and BLOCK_LIST into strings
vcf_df['BLOCK_LIST_2022'] = vcf_df['BLOCK_LIST_2022'].astype(str)
vcf_df['BLOCK_LIST'] = vcf_df['BLOCK_LIST'].astype(str)
tiledb.from_pandas(uri=REPORT_URI,dataframe=vcf_df,mode="ingest",row_start_idx=0)
# In[ ]:


tiledb.__version__


# In[ ]:


import tiledbvcf
tiledbvcf.version


# In[ ]:


import plotly.graph_objects as go


# In[ ]:


metric


# In[ ]:


#just vcf_df
vcf_df_data = metric[metric['name']=='vcf_df'].copy().reset_index(drop=True)


# In[ ]:


fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=vcf_df_data['name'] + '<br>' + vcf_df_data['stage']
    ),
    link=dict(
        source=[0, 1],
        target=[1, 4],
        value=vcf_df_data['value'][[0,1,4]]
    )
)])


# In[ ]:


fig

