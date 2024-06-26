#author: Rumika Mascarenhas
#takes input of the sampleID, affected tissue, unaffected region as the input
#creates a .vcf of variants unique to the affected tissue
#eg: python FilterMutectAndUniq.py 198 198R3 198R5
import sys

#Input files
sampleID = str(sys.argv[1].strip())
affected_sample = str(sys.argv[2].strip())
unaffected_sample = str(sys.argv[3].strip())

#output files
outputfile = affected_sample + "_Mutect_AbsNotInTissueAndGenomic.vcf"
summaryfile = affected_sample + "_Abs_mutect_summary.txt"


print("the affected sample is: " + affected_sample + "\n")
print("the unaffected sample is: " + unaffected_sample + "\n")


#open "database" files
AffectedVCF = open("/work/mtgraovac_lab/Rumika/SeqRound2/MutectRuns/FilteredMutectOut/" + affected_sample + "N_MutectMatchedControl_filtered_PASS.vcf", "r")
UnaffectedVCF = open("/work/mtgraovac_lab/Rumika/SeqRound2/MutectRuns/MutectSingleton/" + unaffected_sample + "N_MutectSingleton_filtered_PASS.vcf", "r")
SEEGR1_TissueVCF = open("/work/mtgraovac_lab/Rumika/SeqRound2/MutectRuns/MutectSingleton/SEEGR2_Tissue_Mutect.vcf", "r")
SEEGR1_GenomicVCF = open("/work/mtgraovac_lab/Rumika/SeqRound2/MutectRuns/MutectSingleton/SEEGR1_Genomic_Mutect.vcf", "r")
SEEGR2_TissueVCF = open("/work/mtgraovac_lab/Rumika/SeqRound2/MutectRuns/MutectSingleton/SEEGR1_" + sampleID + "_Mutect.vcf", "r")
SEEGR2_GenomicVCF = open("/work/mtgraovac_lab/Rumika/SeqRound2/MutectRuns/MutectSingleton/SEEGR2_Genomic_Mutect.vcf", "r")

#1 compare with unaffected.vcf
#2 compare with SEEGR1_Tissue.vcf
#3 compare with SEEGR1_Genomic.vcf
#4 compare with SEEGR2_[sampleID].vcf
#5 compare with SEEGR2_Genomic.vcf

UnaffVarList = []
SEEGR1_Tissue = []
SEEGR1_Genomic = []
SEEGR2_Tissue = []
SEEGR2_Genomic =[]

#make list with variants in unaffected.vcf
for var in UnaffectedVCF:
    if "##" not in var:
        elem = var.split("\t")
        VarLine = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3] + "\t" + elem[4]
        UnaffVarList.append(VarLine.strip())



#make list with variants in SEEGR1_Tissue.vcf
for var in SEEGR1_TissueVCF:
    if "##" not in var:
        elem = var.split("\t")
        VarLine = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3] + "\t" + elem[4]
        SEEGR1_Tissue.append(VarLine.strip())


#make list with variants in SEEGR1_Genomic.vcf
for var in SEEGR1_GenomicVCF:
    if "##" not in var:
        elem = var.split("\t")
        VarLine = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3] + "\t" + elem[4]
        SEEGR1_Genomic.append(VarLine.strip())


#make list with variants in SEEGR2_Tissue.vcf
for var in SEEGR2_TissueVCF:
    if "##" not in var:
        elem = var.split("\t")
        VarLine = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3] + "\t" + elem[4]
        SEEGR2_Tissue.append(VarLine.strip())


#make list with variants in SEEGR2_Genomic.vcf
for var in SEEGR2_GenomicVCF:
    if "##" not in var:
        elem = var.split("\t")
        VarLine = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3] + "\t" + elem[4]
        SEEGR2_Genomic.append(VarLine.strip())



Header = [] #to make sure the format of the new vcf is the same

#Counts

VarCounts = {}

for var2 in AffectedVCF:
    VarCountUnaffected = 0
    VarCountTissue = 0
    VarCountGenomic = 0
    if "##" not in var2:
        elem = var2.split("\t")
        VarLine = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3] + "\t" + elem[4]
        if (VarLine.strip() in SEEGR1_Tissue) or (VarLine.strip() in SEEGR2_Tissue):
            countR1Tissue = SEEGR1_Tissue.count(VarLine.strip())
            countR2Tissue = SEEGR2_Tissue.count(VarLine.strip())
            VarCountTissue = countR1Tissue + countR2Tissue

        if (VarLine.strip() in SEEGR1_Genomic) or (VarLine.strip() in SEEGR2_Genomic):
            countR1Genomic = SEEGR1_Genomic.count(VarLine.strip())
            countR2Genomic = SEEGR2_Genomic.count(VarLine.strip())
            VarCountGenomic = countR1Genomic + countR2Genomic

        if (VarLine.strip() in UnaffVarList):
            VarCountUnaffected = UnaffVarList.count(VarLine.strip())

        #Mutect has 9th column as the unaffected data, which might be causing issue with the annotaiton and filtration
        key = elem[0] + "\t" + elem[1]+ "\t" + elem[2]+ "\t" + elem[3]+ "\t" + elem[4]+ "\t" + elem[5]+ "\t" + elem[6]+ "\t" + elem[7]+ "\t" + elem[8]+ "\t" + elem[9] + "\t" + elem[10]
        VarCounts[key.strip()] = str(VarCountTissue) + "\t" + str(VarCountGenomic) + "\t" + str(VarCountUnaffected) + "\n"
    if "##" in var2:
        #to make sure that the header lines of the vcf are copied
        Header.append(var2.strip())
    if "#CHROM" in var2:
        elem2 = var2.split("\t")
        new_column = elem2[0] + "\t" + elem2[1]+ "\t" + elem2[2]+ "\t" + elem2[3]+ "\t" + elem2[4]+ "\t" + elem2[5]+ "\t" + elem2[6]+ "\t" + elem2[7]+ "\t" + elem2[8]+ "\t" + elem[9] + "\t" + elem2[10]
        Header.append(new_column.strip())


UniqueVar = [] #variants uniq to the affected tissue
CommonVarTissue = [] #variants present in seq round 1 tissue and in seq round 2 tissue
CommonVarGenomic = [] #variants present in seq round 1 genomic and in seq round 2 genomic
CommonVarUnaffected = [] #variants present in unaffected tissue from same sample
CommonVarAll = []

for var3 in VarCounts.keys():
    elem = VarCounts[var3].split("\t")


    CountInTissue = int(elem[0].strip())
    CountInGenomic = int(elem[1].strip())
    CountInUnaffected = int(elem[2].strip())

    #change according to conditions needed
    if (CountInTissue <= 2) and (CountInGenomic == 0) and (CountInUnaffected ==0):
        UniqueVar.append(var3.strip())

    elif (CountInTissue > 0) and (CountInGenomic > 0) and (CountInUnaffected > 0):
        CommonVarAll.append(var3.strip())

    elif (CountInTissue > 0) and (CountInGenomic == 0) and (CountInUnaffected == 0):
        CommonVarTissue.append(var3.strip())

    elif (CountInTissue == 0) and (CountInGenomic > 0) and (CountInUnaffected == 0):
        CommonVarGenomic.append(var3.strip())

    elif (CountInTissue == 0) and (CountInGenomic == 0) and (CountInUnaffected > 0):
        CommonVarUnaffected.append(var3.strip())



#testing
#print(len(UniqueVar))
#print(len(CommonVarGenomic))
#print(len(CommonVarTissue))
#print(len(CommonVarUnaffected))
#print(len(CommonVarAll))



#Output 
output = open(outputfile, "w")
outputsummary = open(summaryfile , "w")


for i in range(len(Header)):
    output.write(Header[i] + "\n")

for i in range(len(UniqueVar)):
    output.write(UniqueVar[i] + "\n")


#Summary file
outputsummary.write("Number of variants present in Genomic samples:\t" + str(len(CommonVarGenomic)) +"\n")
outputsummary.write("Number of variants present in Tissue samples:\t" + str(len(CommonVarTissue)) +"\n")
outputsummary.write("Number of variants present in Unaffected sample:\t" + str(len(CommonVarUnaffected)) +"\n")
outputsummary.write("Number of variants unique to Affected sample:\t" + str(len(UniqueVar)) +"\n")

