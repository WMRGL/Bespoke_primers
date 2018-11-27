import os
import csv
import autoprimer.autoprimer as ap
import requests
import sys
from pybedtools import BedTool


class Variant:
    def __init__(self, chromosome, start, build, **kwargs):
        self.filename = ""
        self.gene = ""
        self.strand = ""
        self.build = build
        self.chromosome = chromosome
        self.start = start
        self.end = ""
        self.lenght = ""
        self.ref = ""
        self.alt = ""
        self.inheritance = ""
        self.condition = ""
        self.hgvsc = ""
        self.hgvsp = ""
        self.zygosity = ""
        self.pathogenicity = ""
        self.contribution = ""
        self.depth = ""
        self.af_max = ""
        for key, value in kwargs.items():
            setattr(self, key, value)


def get_surrounding_sequence(variant):
    start = int(variant.start) - 501
    stop = int(variant.start) + 500
    variant_bed = BedTool(variant.chromosome + " " + str(start) + " " + str(stop), from_string=True)
    # set the genome to use
    if variant.build == 'GRCh37':
        ucsc_fasta = BedTool("/media/sf_S_DRIVE/genomic_resources/primer_design/hg19.fa")
        #ucsc_fasta = BedTool("/srv/primer_design/s_drive/hg19.fa")
    elif variant.build == 'GRCh38':
        ucsc_fasta = BedTool("/media/sf_S_DRIVE/genomic_resources/primer_design/hg38.fa")
        #ucsc_fasta = BedTool("/srv/primer_design/s_drive/hg38.fa")
    # use pybedtools API to return the sequence
    genomic_region = variant_bed.sequence(fi=ucsc_fasta, tab=True)
    bedtools_result = open(genomic_region.seqfn).read()
    raw_sequence = bedtools_result.strip().split('\t')[1]
    return raw_sequence.upper()


def look_up_strand(gene_symbol):
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens/" + gene_symbol + "?expand=1"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    if repr(decoded['strand']) == "1":
        strand = "Forward (+)"
    elif repr(decoded['strand']) == "-1":
        strand = "Reverse (-)"
    return strand


def design_from_coord(variant, options):
    print('chromosome = ' + variant.chromosome)
    print('coord = ' + str(variant.start))
    sequence = get_surrounding_sequence(variant)
    input_sequence = ap.InputSequence(variant.hgvsc, sequence, gene_name=variant.gene, chrom_number=variant.chromosome,
                                      genomic_coords=(int(variant.start)-500, int(variant.start)+500), strand="+")
    target = ap.TargetRegion(variant.chromosome + ":" + str(variant.start), sequence, int(variant.start)-500, int(variant.start)+500, 450, input_sequence)
    print(type(target))
    input_sequence.set_snps_bed(variant.build)
    input_sequence.target_regions.append(target)

    target_regions = input_sequence.target_regions
    for target in target_regions:
        target.set_snps()
        target.mask_sequence(options['max_avhet'])
        print(target.masked_sequence)
        target.set_primers(options['min_product_size'], options['max_product_size'], options['primer_opt_size'],
                           options['primer_min_size'], options['primer_max_size'], options['primer_opt_tm'],
                           options['primer_min_tm'], options['primer_max_tm'], options['primer_min_gc'],
                           options['primer_max_gc'])
        primers = target.primers
        for primer in primers:
            primer.set_snps(target)
            print(primer.forward_seq, primer.forward_genomic_coords, primer.reverse_seq, primer.reverse_genomic_coords)
    return input_sequence



def write_to_csv(input_sequence, variant):
    """
    Writes primers that have been designed to a CSV file
    """
    # inputSequence ID used as the filename
    filename = variant.chromosome + '-' + str(variant.start) + '.csv'
    targets = input_sequence.target_regions
    #with open('/srv/primer_design/s_drive/designs/' + filename, 'w') as csvfile:
    with open('/media/sf_S_DRIVE/genomic_resources/primer_design/designs/' + filename, 'w') as csvfile:
        f = csv.writer(csvfile, delimiter=',',
                       quotechar=',', quoting=csv.QUOTE_MINIMAL)
        f.writerow(['Gene', 'Strand', 'Target', 'Product size', 'Forward primer sequence', 'Genomic Coords', 'Forward TM',
                    'Forward GC %', 'Forward SNPs', 'Reverse primer sequence', 'Genomic Coords', 'Reverse TM',
                    'Reverse GC %', 'Reverse SNPs'])
        for target in targets:
            primer_list = target.primers
            # Primer temperatures and GC% rounded to 2 decimal places
            for primer in primer_list:
                forward_snps = ''
                reverse_snps = ''
                for snp in primer.forward_snps:
                    forward_snps = forward_snps + snp.snp_id + ' (' + str(round(snp.av_het, 4)) + ') '
                for snp in primer.reverse_snps:
                    reverse_snps = reverse_snps + snp.snp_id + ' (' + str(round(snp.av_het, 4)) + ') '
                f.writerow([input_sequence.gene_name, variant.strand, target.target_id, primer.product_size, primer.forward_seq,
                            input_sequence.chrom_number + ":" + str(primer.forward_genomic_coords[0]) + "-" + str(
                                primer.forward_genomic_coords[1]), round(primer.forward_tm, 2),
                            round(primer.forward_gc, 2), forward_snps, primer.reverse_seq,
                            input_sequence.chrom_number + ":" + str(primer.reverse_genomic_coords[0]) + "-" + str(
                                primer.reverse_genomic_coords[1]),
                            round(primer.reverse_tm, 2), round(primer.reverse_gc, 2), reverse_snps])

def write_to_bed(input_sequence, variant):
    """
    Writes primers that have been designed to a CSV file
    """
    # inputSequence ID used as the filename
    filename = variant.chromosome + '-' + str(variant.start) + '.bed'
    targets = input_sequence.target_regions
    #with open('/srv/primer_design/s_drive/designs/' + filename, 'w') as csvfile:
    with open('/media/sf_S_DRIVE/genomic_resources/primer_design/designs/' + filename, 'w') as csvfile:
        f = csv.writer(csvfile, delimiter='\t',
                       quotechar=';', quoting=csv.QUOTE_MINIMAL)
        f.writerow(['track name="' + filename + '" description=' + '"Primers designed for' + filename +
                    '" visibility=2 itemRgb="On"'])
        for target in targets:
            f.writerow([input_sequence.chrom_number, target.seq_start + target.overhang - 30,
                        target.seq_stop - target.overhang + 30, target.target_id, 0, input_sequence.strand,
                        target.seq_start + target.overhang - 30, target.seq_stop - target.overhang + 30,
                        '255,0,0'])
        for target in targets:
            for primer in target.primers:
                f.writerow([input_sequence.chrom_number, primer.forward_genomic_coords[0] - 1,
                            primer.reverse_genomic_coords[1], target.target_id, 0, input_sequence.strand,
                            primer.forward_genomic_coords[1], primer.reverse_genomic_coords[0] - 1, '0,0,255'])


def design_by_coordinate(genome_build, chromosome, coordinate, max_avhet, min_product_size, max_product_size, primer_opt_size,
                             primer_min_size, primer_max_size, primer_opt_tm, primer_min_tm, primer_max_tm,
                             primer_min_gc, primer_max_gc):
    variant = Variant(chromosome, coordinate, genome_build)
    options = {'max_avhet': float(max_avhet),
               'min_product_size': int(min_product_size),
               'max_product_size': int(max_product_size),
               'primer_opt_size': int(primer_opt_size),
               'primer_min_size': int(primer_min_size),
               'primer_max_size': int(primer_max_size),
               'primer_opt_tm': float(primer_opt_tm),
               'primer_min_tm': float(primer_min_tm),
               'primer_max_tm': float(primer_max_tm),
               'primer_min_gc': float(primer_min_gc),
               'primer_max_gc': float(primer_max_gc)}
    input_sequence = design_from_coord(variant, options)
    write_to_csv(input_sequence, variant)
    write_to_bed(input_sequence, variant)
    return input_sequence