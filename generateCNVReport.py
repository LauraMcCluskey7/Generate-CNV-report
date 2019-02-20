import pandas
import numpy
import allel
import argparse

dataFrame=pandas.DataFrame()

parser=argparse.ArgumentParser()
parser.add_argument('--runid')
parser.add_argument('--output')
args=parser.parse_args()
runid=args.runid
output_file=args.output

samples_list=[]
Exome_depth_metrics= pandas.read_table("/data/results/"+ runid +"/IlluminaTruSightCancer/" + runid+ "_ExomeDepth_Metrics.txt") 
num_rows_exome_depth_metrics= Exome_depth_metrics.shape[0]

#extract sample names from exome_depth_metrics file and create sample list
row=0
while (row<num_rows_exome_depth_metrics):
    samples_bam= Exome_depth_metrics.iloc[row,0].split("_")
    sample= samples_bam[7].split(".")
    sample_2=sample[0]
    samples_list.append(sample_2)
    row=row+1

#Loop through all samples and read in the relevant files
for sample in samples_list: 
    bed_file= pandas.read_table("/data/results/" +runid+ "/IlluminaTruSightCancer/IlluminaTruSightCancer_CustomROI_b37.bed")
    depth_of_coverage_summary=pandas.read_table("/data/results/" + runid+ "/IlluminaTruSightCancer/" +sample+ "/"+ runid+ "_"+ sample+ "_DepthOfCoverage.sample_summary")
    depth_of_coverage_summary_mean=depth_of_coverage_summary['mean'][1]
    
    #Create dataframe from manta vcf
    manta_dataframe=allel.vcf_to_dataframe("/data/results/" +runid+ "/IlluminaTruSightCancer/" +sample +"/" +runid + "_" + sample + "_sv_filtered.vcf.gz" , fields=['*'])
    if ((manta_dataframe)is not None):
        manta_dataframe_2=manta_dataframe.iloc[:,[0,1,11,3,9]]
        manta_dataframe_2.columns=['CHROM', 'POS', 'END', 'REF', 'ALT_1']
        manta_dataframe_2['Regions']= "-"
        manta_dataframe_2["QC"]="PASS"
        manta_dataframe_2['method']="Manta"
    else:
        data= [{'CHROM':'NA', 'POS':'NA', 'END':'NA', 'REF':'NA', 'ALT_1':'NA', 'Regions':'NA'}]
        manta_dataframe_2=pandas.DataFrame(data)
        manta_dataframe_2=manta_dataframe_2.iloc[:,[1,3,2,4,0,5]]
        manta_dataframe_2["QC"]="PASS"
        manta_dataframe_2['method']="Manta"

    #Create dataframe from Exome depth vcf
    Exome_dataFrame=allel.vcf_to_dataframe("/data/results/" + runid+ "/IlluminaTruSightCancer/" +sample+ "/" +runid+ "_" + sample + "_cnv.vcf.gz", fields=['*'])
    if ((Exome_dataFrame)is not None):
        Exome_dataFrame_2=Exome_dataFrame.iloc[:,[0,1,8,3,4,10]]
        Exome_dataFrame_2["QC"]="PASS"
        Exome_dataFrame_2['method']="Exome"
    else:
        Exome_data= [{'CHROM':'NA', 'POS':'NA', 'END':'NA', 'REF':'NA', 'ALT_1':'NA', 'Regions':'NA'}]
        Exome_dataFrame_2=pandas.DataFrame(Exome_data)
        Exome_dataFrame_2=Exome_dataFrame_2.iloc[:,[1,3,2,4,0,5]]
        Exome_dataFrame_2["QC"]="PASS"
        Exome_dataFrame_2['method']="Exome"
   
    #Add columns for depth and sampleid
    if (((Exome_dataFrame_2)is not None) and ((manta_dataframe_2)is not None)):
        Manta_Exome= Exome_dataFrame_2.append(manta_dataframe_2)
        Manta_Exome['depth']= depth_of_coverage_summary_mean
        Manta_Exome['sampleid']= sample
    elif (((Exome_dataFrame_2)is None) and ((manta_dataframe_2)is not None)):
        Manta_Exome= manta_dataframe_2
        Manta_Exome['depth']= depth_of_coverage_summary_mean
        Manta_Exome['sampleid']= sample
    elif (((Exome_dataFrame_2)is not None) and ((manta_dataframe_2)is None)):
        Manta_Exome= Exome_dataFrame_2
        Manta_Exome['depth']= depth_of_coverage_summary_mean
        Manta_Exome['sampleid']= sample
    else:
        Manta_Exome= None

    #If end position is -1, change value to start_position+1   
    if (Manta_Exome is not None):
        row_num=0
        Manta_exome_num_rows=Manta_Exome.shape[0]
        while (row_num<Manta_exome_num_rows):
            if (Manta_Exome.iloc[row_num,2]!="NA"):
                Manta_Exome.iloc[row_num,2]= int(Manta_Exome.iloc[row_num,2])
                if (Manta_Exome.iloc[row_num, 2]== (-1)):
                    Manta_Exome.iloc[row_num, 2]= Manta_Exome.iloc[row_num,1]+1
            row_num=row_num +1 

    #Change value in QC column if depth<160 or R2<0.98 for sample
    if (Manta_Exome is not None):    
    
        row=0 
        num_rows_exome_depth_metrics= Exome_depth_metrics.shape[0]

        while (row<num_rows_exome_depth_metrics):
            samples_bam= Exome_depth_metrics.iloc[row,0].split("_")
            sample_2= samples_bam[7].split(".")
            sample_3=sample_2[0]
            if (sample_3 ==sample):
                Manta_Exome['QC']= Manta_Exome['QC'].astype(str)
                a=0
                for i in Manta_Exome['depth']:
                    print(Manta_Exome.iloc[a,6])
                    Manta_Exome.iloc[a,6]
                    if (i<160 and Exome_depth_metrics.iloc[row, 2]<0.98):
                        Manta_Exome.iloc[a,6]="R2<0.98; Depth<160"
                    elif (i<160):
                        Manta_Exome.iloc[a,6]="Depth<160"
                    a=a+1
            row=row+1

    
    # Add the Gene column by comparing chrom, start position and end position to those in bed file
    if (Manta_Exome is not None):
        num_rows_Manta_Exome= Manta_Exome.shape[0]
        num_rows_bed_file= bed_file.shape[0]

        a=0
        b=0

        gene_list=[]

        while (a<num_rows_Manta_Exome):
            b=0
            while (b<num_rows_bed_file):
                if (type(bed_file.iloc[b,0])!= str):
                    bed_file.iloc[b,0]=int(bed_file.iloc[b,0])
                else:
                    Manta_Exome.iloc[a,0]=str(Manta_Exome.iloc[a,0])
                bed_file.iloc[b,1]=int(bed_file.iloc[b,1])
                bed_file.iloc[b,2]=int(bed_file.iloc[b,2])
                if(type(Manta_Exome.iloc[a,0]) != str):
                        Manta_Exome.iloc[a,0]= int(Manta_Exome.iloc[a,0])
                else:
                    bed_file.iloc[b,0]=str(bed_file.iloc[b,0])
                if(Manta_Exome.iloc[a,1] != "NA"):
                    Manta_Exome.iloc[a,1]= int(Manta_Exome.iloc[a,1])
                if(Manta_Exome.iloc[a,2] != "NA"):    
                    Manta_Exome.iloc[a,2]= int(Manta_Exome.iloc[a,2])
                if (((Manta_Exome.iloc[a,0]== bed_file.iloc[b,0]) and(Manta_Exome.iloc[a,2]>= bed_file.iloc[b,2]) and (Manta_Exome.iloc[a,1]<= bed_file.iloc[b,1])) or ((Manta_Exome.iloc[a,0]== bed_file.iloc[b,0]) and(Manta_Exome.iloc[a,1]>= bed_file.iloc[b,1]) and (Manta_Exome.iloc[a,2]<= bed_file.iloc[b,2])) or ((Manta_Exome.iloc[a,0]== bed_file.iloc[b,0]) and (Manta_Exome.iloc[a,1]<= bed_file.iloc[b,1]) and (Manta_Exome.iloc[a,2]<= bed_file.iloc[b,2])and (Manta_Exome.iloc[a,2] >= bed_file.iloc[b,1])) or ((Manta_Exome.iloc[a,0]== bed_file.iloc[b,0]) and (Manta_Exome.iloc[a,1]<= bed_file.iloc[b,2]) and (Manta_Exome.iloc[a,1]>= bed_file.iloc[b,1]) and (Manta_Exome.iloc[a,2]>=bed_file[b,2]))) :
                    bed_file.iloc[b,3]=str(bed_file.iloc[b,3])
                    gene= bed_file.iloc[b,3].split(".")
                    if (len(gene_list)<=a):
                        gene_list.append(gene[0])
                
            
                b=b+1
            if (len(gene_list)<=a):
                gene_list.append("NA")
    
            a=a+1

    if (Manta_Exome is not None):
        Manta_Exome['gene']= gene_list

    #Reorder the columns and change column names
    if (Manta_Exome is not None):
        order= [9,7,0,1,2,3,4,5,6,8,10]
        Manta_Exome=Manta_Exome[Manta_Exome.columns[order]]
        Manta_Exome.columns=['sample', 'method','chr', 'start', 'end', 'ref', 'type', 'regions', 'qc', 'depth', 'gene']

    #Add sample dataframe to main dataframe containing all samples
    if (Manta_Exome is not None):
        dataFrame=dataFrame.append(Manta_Exome)
        
#output the dataframe containing the CNVs for all sample as a CSV
dataFrame.to_csv(output_file, index=False)


