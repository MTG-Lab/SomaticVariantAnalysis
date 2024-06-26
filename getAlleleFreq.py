#getAlleleFreq.py

import sys

#takes this from the command line
Affected = open(sys.argv[2], "r")
Unaffected = open(sys.argv[3], "r")
Genomic = open(sys.argv[1], "r")

a = sys.argv[2].split(".")
filename = a[0].strip() + "AF_out.txt"

VarGenomic = {}
VarAffected = {}
VarUnaffected = {}

for var in Genomic:
    if "CHROM" not in var:
        elem = var.split("\t")
        AD = elem[5].split(",")
        TotalDepth = elem[4].strip()
        AlleleDepth = AD[1].strip()
        AF = float(AlleleDepth)/float(TotalDepth)

        VarID = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3]
        VarGenomic[VarID] = str(AF)

for var in Affected:
    if "CHROM" not in var:
        elem = var.split("\t")
        AD = elem[5].split(",")
        TotalDepth = elem[4].strip()
        AlleleDepth = AD[1].strip()
        AF = float(AlleleDepth)/float(TotalDepth)

        VarID = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3]
        VarAffected[VarID] = str(AF)


for var in Unaffected:
    if "CHROM" not in var:
        elem = var.split("\t")
        AD = elem[5].split(",")
        TotalDepth = elem[4].strip()
        AlleleDepth = AD[1].strip()
        AF = float(AlleleDepth)/float(TotalDepth)

        VarID = elem[0] + "\t" + elem[1] + "\t" + elem[2] + "\t" + elem[3]
        VarUnaffected[VarID] = str(AF)


output = open(filename, "w")

output.write("CHROM\t" + "POS\t" + "REF\t" + "ALT\t" + "GenomicAF\t" + "AffectedAF\t" + "UnaffectedAF\n")

for x in VarGenomic.keys():
    if (x in VarAffected.keys()) and (x in VarUnaffected.keys()):
        output.write(x + "\t")
        output.write(VarGenomic[x] + "\t")
        output.write(VarAffected[x] + "\t")
        output.write(VarUnaffected[x] + "\n")

