#!/usr/bin/env python3

from collections import defaultdict
import numpy as np
RTE_tab = r'/blue/soltis/shan158538/Methylation/OutPut/TE-gene_distance/repeat_annotation_combined_LINE.final_4.gff' #r'D:\Programming\R_working\Tyrran_paper\TEclassification7.tab'
gff3 = r'/blue/soltis/shan158538/Methylation/OutPut/TE-gene_distance/Tdub.V1.rm.RENAME.gff' #r'sbs_20141217-Ha412v1r1.gff3'

# I changed the code of how to extract the te_id and the chromosome name.
def getRTEmiddleCoordinates(RTE_tab):
    RTE_middle_point = defaultdict(list) #chromosome:[middle points of RTEs...]
    middle_point_TE_id = {} # chromosome:{middel point:TE_id}
    with open(RTE_tab) as tab:
        for i, lines in enumerate(tab):
            if i != 0:
                sp = lines.split("\t")
                te_id, chromosome, middle = sp[0].split(";")[0], \
                                            sp[1], \
                                            int((int(sp[2]) + int(sp[3]))/2)
                RTE_middle_point[chromosome].append(middle)

                if chromosome not in middle_point_TE_id:
                    middle_point_TE_id[chromosome] = {}

                middle_point_TE_id[chromosome][middle] = te_id
    return [RTE_middle_point, middle_point_TE_id]

# I changed the minimum length of a gene from 1,000 to 0.
def getGenePositionFromGff3(gff3, min_gene_len = 0):
    Gene_middle_point = defaultdict(list)  # chromosome:[middle points of Gene...]
    middle_point_Gene_id = {}  # chromosome:{middel point:Gene id}
    cnt = 0
    with open(gff3) as gff:
        for line in gff:
            if not line.startswith("##"):
                sp = line.split("\t")
                if sp[2] == "gene" and abs(int(sp[4]) - int(sp[3])) >= min_gene_len:
                    gene_id, chromosome, middle = sp[-1].split(";")[0], \
                                                sp[0], \
                                                int((int(sp[4]) + int(sp[3])) / 2)
                    Gene_middle_point[chromosome].append(middle)

                    if chromosome not in middle_point_Gene_id:
                        middle_point_Gene_id[chromosome] = {}

                    middle_point_Gene_id[chromosome][middle] = [gene_id, min(int(sp[3]), int(sp[4])), max(int(sp[3]), int(sp[4]))]
                    cnt += 1
    print("Number of genes:", cnt)
    return [Gene_middle_point, middle_point_Gene_id]

def getClosestGene(RTE_middle_point, Gene_middle_point):
    rte_gene_middle_points = defaultdict(list) #chromosome:[[rte middle, closest gene middle],[]]
    for chromosomes in Gene_middle_point:
        print("Chromosome {0} of {1}".format(chromosomes, len(Gene_middle_point)))
        sorted_genes_chromosomes = sorted(Gene_middle_point[chromosomes])
        for rte_middle_point in RTE_middle_point[chromosomes]:
            idx = np.searchsorted(sorted_genes_chromosomes, rte_middle_point)
            if idx == 0: # inserted in the beginning of the list
                distance = abs(sorted_genes_chromosomes[0] - rte_middle_point)
                closest_gene_middle = sorted_genes_chromosomes[0]
            elif idx == len(sorted_genes_chromosomes): # inserted in the end of the list (then index will be == len(sorted_genes_chromosomes)
                distance = abs(sorted_genes_chromosomes[idx - 1] - rte_middle_point)
                closest_gene_middle = sorted_genes_chromosomes[idx - 1]
            else:
                distance = min( abs(sorted_genes_chromosomes[idx - 1] - rte_middle_point),
                                    abs(sorted_genes_chromosomes[idx] - rte_middle_point)
                                    )
                ind = [abs(sorted_genes_chromosomes[idx - 1] - rte_middle_point),
                                    abs(sorted_genes_chromosomes[idx] - rte_middle_point)].index(distance)

                if ind == 0:
                    closest_gene_middle = sorted_genes_chromosomes[idx - 1]
                elif ind == 1:
                    closest_gene_middle = sorted_genes_chromosomes[idx]

            rte_gene_middle_points[chromosomes].append([rte_middle_point, closest_gene_middle, distance])
    return rte_gene_middle_points
def main():
    RTE_middle_point, middle_point_TE_id = getRTEmiddleCoordinates(RTE_tab)
    Gene_middle_point, middle_point_Gene_id = getGenePositionFromGff3(gff3,min_gene_len=0)
    closest_genes = getClosestGene(RTE_middle_point, Gene_middle_point)

    with open("closest_gene.tab", "w") as out:
        out.write("\t".join(["uniqueTE_id", "TE.middle.point", "Gene.id", "Gene.start", "Gene.end", "Gene.middle.point", "Distance"]) + "\n")
        for chromosomes in closest_genes:
            #print(closest_genes[chromosomes])
            #print(middle_point_TE_id[chromosomes])
            for te_middle, gen_midle, distance in closest_genes[chromosomes]:
                te_id, gene_id = middle_point_TE_id[chromosomes][te_middle], middle_point_Gene_id[chromosomes][gen_midle][0]
                gene_start, gene_end = middle_point_Gene_id[chromosomes][gen_midle][1], middle_point_Gene_id[chromosomes][gen_midle][2]
                out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(te_id, te_middle, gene_id, gene_start, gene_end, gen_midle, distance))

main()

#print(sorted([1,4,6,2,5]))
#print(np.searchsorted(sorted([1,4,6,2,5]), 3))