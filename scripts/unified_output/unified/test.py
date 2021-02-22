import csv

psi_file = "/home/alex/Documents/Studium_Bioinformatik/Master/Semester3/FortgeschrittenenPrakitkum/unifyASOutputs/data/majiq-output/two_event_50M_psi_voila.tsv"
voila_file = "/home/alex/Documents/Studium_Bioinformatik/Master/Semester3/FortgeschrittenenPrakitkum/unifyASOutputs/data/majiq-output/two_event_50M_voila.tsv"
merged_file = "/home/alex/Documents/Studium_Bioinformatik/Master/Semester3/FortgeschrittenenPrakitkum/unifyASOutputs/data/majiq-output/two_event_50M_merged_voila.tsv"

header = ["gene_id","lsv_id","num_junctions","num_exons","de_novo_junctions","seqid","strand","junctions_coords","exons_coords","ir_coords","A3SS","A5SS","ES"]


psi_dict={}
with open(psi_file, "r") as f:
    # read header
    f.readline()
    for line in f:
        line = line.strip().split("\t")
        idx = line[1]
        a5 = line[5]
        a3 = line[6]
        es = line[7]
        psi_dict[idx] = [a3, a5, es]
f.close()

voila_dict = {}
with open(voila_file, 'r') as f:
    while True:
        line = f.readline()
        if not line.startswith("#"): break
    h = line.split("\t")
    for line in f:
        line=line.strip().split("\t")
        idx = line[1]
        voila_dict[idx] = [line[0], line[1], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12]]
#gene_id	lsv_id	mean_psi_per_lsv_junction	stdev_psi_per_lsv_junction	lsv_type	num_junctions	num_exons	de_novo_junctions	seqid	strand	junctions_coords	exons_coords	ir_coords	ucsc_lsv_link
#0           1                                                                                   5             6             7                8       9          10               11            12
with open(merged_file, 'w') as out:
    tsv_writer = csv.writer(out, delimiter='\t')
    tsv_writer.writerow(header)
    for idx, values in voila_dict.items():
        values.extend(psi_dict[idx])
        tsv_writer.writerow(values)




