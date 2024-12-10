#Write a Python script named pick_fa.py to extract nucleic acid sequences from genbank files, such as "SRR5935758_all.1.bgc.gbk", "SRR5935758_all.2.bgc.gbk", "SRR5946523_all.1.bgc.gbk"... and other genbank files.
from Bio import SeqIO
gbk_filename = "SRR5935758_all.1.bgc.gbk"
faa_filename = "SRR5935758_all.1.bgc.fna"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    output_handle.write(">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))

output_handle.close()
input_handle.close()
gbk_filename = "SRR5935758_all.2.bgc.gbk"
faa_filename = "SRR5935758_all.2.bgc.fna"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    output_handle.write(">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))

output_handle.close()
input_handle.close()
gbk_filename = "SRR5946523_all.1.bgc.gbk"
faa_filename = "SRR5946523_all.1.bgc.fna"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    output_handle.write(">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))

output_handle.close()
input_handle.close()