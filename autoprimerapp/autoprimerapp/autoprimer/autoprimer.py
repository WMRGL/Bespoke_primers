#! /usr/bin/python


# To do:
#     Add in a check to make sure the nucleotide being masked matches one of the nucleotides in the SNP.
#     Add in functionality to adjust primer3 options
#     Accept input as gene name?
#     Get isPCR to work and incorporate into project.
#     Determine fields used in Sarah's primer database so correct output can be given.
#     Error handling.
#     User interface? Framework?
#     Testing

import sys
import csv
import os
import primer3
import math
import re
from Bio import Entrez, SeqIO, SearchIO
from Bio.Seq import Seq
from pybedtools import BedTool


class InputSequence:
    """
    InputSequence instances initalised with the ID and sequence extracted from NCBI using the supplied NG number.
    Other instance variables record information about genomic location and a list of all SNPs (version 147) within
    the region in bed format. This list is used to mask common SNPs when desgining primers and for identifying other
    primer site SNPs.
    A list of all exons in the gene are stored and from these a list of target regions are determined.
    """

    def __init__(self, refseq_id, sequence, **kwargs):
        self.refseq_id = refseq_id
        self.sequence = sequence
        self.sequence_length = len(sequence)
        self.gene_name = ""
        self.chrom_number = None
        self.genomic_coords = None
        self.strand = None
        self.snps_bed = None
        self.exons = []
        self.target_regions = []
        for key, value in kwargs.items():
            setattr(self, key, value)

    @staticmethod
    def blat_search(fasta_file, genome_build):
        """
        Performs a BLAT search of the input sequence.
        Please note, this uses a client search of a local server. The local server needs to be
        initalised by running blatServer.py prior to running this script. The hits are output
        to a .psl file which is saved in the output folder.
        """
        # the blat search requires an input file, so the sequence is first saved to disk in FASTA format.
        f = open("input/fasta.fa", "w")
        f.write(fasta_file)
        f.close()
        # used minIdentity and minScore parameters to reduce the number of positive hits.
        if genome_build == 'GRCh37':
            os.system('gfClient localhost 17777 / -minIdentity=100 -minScore=100' + ' ' +
                      "input/fasta.fa" + " output/gfOutput.psl")
            # os.system('blat ' + '/opt/genome/hg19.2bit' + ' ' + 'input/fasta.fa' + " output/gfOutput.psl")
        elif genome_build == 'GRCh38':
            os.system('gfClient localhost 17778 / -minIdentity=100 -minScore=100' + ' ' +
                      "input/fasta.fa" + " output/gfOutput.psl")

    def set_genomic_location(self, psl_file_path):
        """
        Sets the chromosomeNumber, genomicCoords and strand by parsing the .psl BLAT output.
        Uses SearchIO from the Biopython library to parse the .psl file.
        NB: SearchIO is an experimental submodule which may undergo significant changes prior to official release
        """
        coords_set = False
        qresult = SearchIO.read(psl_file_path, 'blat-psl')
        for hit in qresult:
            for hsp in hit:
                # checks that the hit is a complete match with the given sequence. The maximum size for BLAT searches
                # is ~ 45000 bases, so if the sequence length is great than this the match size will not equal the
                # sequence length. Therefore hsp.match_num > 40000 is used to catch any very large sequence matches.
                # Tested on DDX41 NG_046846 and no exact blat match found due to 1 base being wrong.. So added in
                # hsp.match_num > qresult.seq_len * 0.99 to catch any very close but not exact matches.
                if hsp.match_num == qresult.seq_len or hsp.match_num > 40000 or hsp.match_num > qresult.seq_len * 0.99:
                    self.chrom_number = hit.id
                    # Used +1 to capture location of first base - checked start and stop coords on UCSC genome browser.
                    self.genomic_coords = (hsp.hit_start + 1, hsp.hit_start + self.sequence_length)
                    for hspFragment in hsp:
                        self.strand = hspFragment.query_strand
                        # Strand recorded as '+' or '-'
                        if self.strand == 1:
                            self.strand = "+"
                        else:
                            self.strand = "-"
                        coords_set = True
        if not coords_set:
            raise ValueError("No exact sequence match found from BLAT search")

    def set_exons(self, refseq):
        """
        Sets the each exon from the reference sequence. The local start and stop positions are the locations of the
        exons within the reference sequence. The genomic locations for the exons are calculated using these local
        positions and the genomic location of the reference sequence as determined from the BLAT search.
        """
        exons = []
        for feature in refseq.features:
            if feature.type == "exon":
                self.gene_name = feature.qualifiers.get('gene')[0]
                number = feature.qualifiers.get('number')[0]
                if self.strand == "+":
                    loc_start = feature.location._start - 1
                    loc_stop = feature.location._end - 1
                    genom_start = loc_start + self.genomic_coords[0]
                    genom_stop = loc_stop + self.genomic_coords[0]
                else:
                    loc_start = feature.location._start + 1
                    loc_stop = feature.location._end + 1
                    genom_start = self.genomic_coords[1] - loc_start
                    genom_stop = self.genomic_coords[1] - loc_stop
                exon = Exon(number, loc_start, loc_stop, genom_start, genom_stop)
                exons.append(exon)
        exons[-1].final_exon = True
        self.exons = exons

    # function may need refactoring...
    def set_target_regions(self, max_size, overhang):
        """
        Determines the target regions from the exons. If an exon is over the maximum size then this will be split into
        smaller target regions and if multiple small exons are close together then these will be grouped as a single
        target providing the target size remains less than the maximum size. The overhang represents the sequence
        either side of the target to which the primers can be designed.
        """
        targets = []
        for i, exon in enumerate(self.exons):
            # if start position of stop position of exon - start postion of exon > max_size, split into smaller chucks
            # maxSize currently hard coded at 490 but this could be determined by the user
            if exon.length > max_size:
                exon.processed = True
                number_of_fragments = math.ceil(exon.length / max_size)  # rounded up to nearest 1
                size_of_each_fragment = math.ceil(exon.length / number_of_fragments)  # rounded up to nearest 1
                loc_start = exon.loc_start
                loc_stop = exon.loc_start
                genom_start = exon.genom_start
                genom_stop = exon.genom_start
                count = 1
                for _ in range(number_of_fragments):
                    loc_stop += size_of_each_fragment
                    # creates a cropped version of the input sequence around the target region
                    cropped_seq = self.sequence[loc_start - overhang:loc_stop + overhang]
                    if self.strand == "+":
                        genom_stop += size_of_each_fragment
                        cropped_seq_genom_start = genom_start - overhang
                        cropped_seq_genom_stop = genom_stop + overhang
                    else:
                        genom_stop -= size_of_each_fragment
                        cropped_seq_genom_start = genom_start + overhang
                        cropped_seq_genom_stop = genom_stop - overhang
                    target_id = str(exon.number) + "." + str(count)
                    target_region = TargetRegion(target_id, cropped_seq, cropped_seq_genom_start,
                                                 cropped_seq_genom_stop, overhang, self)
                    targets.append(target_region)
                    loc_start = loc_stop
                    genom_start = genom_stop
                    count += 1

            # if start position of exon + stop position of next exon <= max_size then group exons together
            elif not exon.final_exon:
                if (self.exons[i + 1].loc_stop - exon.loc_start) <= max_size and not exon.processed:
                    exon.processed = True
                    self.exons[i + 1].processed = True
                    stop_exon = self.exons[i + 1]
                    loc_start = exon.loc_start
                    loc_stop = self.exons[i + 1].loc_stop
                    count = 2
                    # while loop checks to see if next exons should also be included in same target region
                    if not self.exons[i + 1].final_exon:
                        while (self.exons[i + count].loc_stop - loc_start) <= max_size and count <= (
                                    len(self.exons) - int(re.findall("\d+", exon.number)[0])):
                            self.exons[i + count].processed = True
                            loc_stop = self.exons[i + count].loc_stop
                            stop_exon = self.exons[i + count]
                            if self.exons[i + count].final_exon:
                                break
                            count += 1
                    target_id = str(exon.number) + "-" + str(stop_exon.number)
                    cropped_seq = self.sequence[loc_start - overhang:loc_stop + overhang]
                    if self.strand == "+":
                        genom_start = exon.genom_start
                        genom_stop = stop_exon.genom_stop
                        cropped_seq_genom_start = genom_start - overhang
                        cropped_seq_genom_stop = genom_stop + overhang
                    else:
                        genom_start = exon.genom_start
                        genom_stop = stop_exon.genom_stop
                        cropped_seq_genom_start = genom_start + overhang
                        cropped_seq_genom_stop = genom_stop - overhang
                    target_region = TargetRegion(target_id, cropped_seq, cropped_seq_genom_start,
                                                 cropped_seq_genom_stop, overhang, self)
                    targets.append(target_region)

            # if the exon is not too large or not close enough to another small exon then the exon is treated as a
            # target region on its own.
            if not exon.processed:
                target_id = exon.number
                loc_start = exon.loc_start
                loc_stop = exon.loc_stop
                cropped_seq = self.sequence[loc_start - overhang:loc_stop + overhang]
                if self.strand == "+":
                    genom_start = exon.genom_start
                    genom_stop = exon.genom_stop
                    cropped_seq_genom_start = genom_start - overhang
                    cropped_seq_genom_stop = genom_stop + overhang
                else:
                    genom_start = exon.genom_start
                    genom_stop = exon.genom_stop
                    cropped_seq_genom_start = genom_start + overhang
                    cropped_seq_genom_stop = genom_stop - overhang
                target_region = TargetRegion(target_id, cropped_seq, cropped_seq_genom_start, cropped_seq_genom_stop,
                                             overhang, self)
                targets.append(target_region)

        self.target_regions = targets

    def set_snps_bed(self, genome_build):
        """
        Uses the genomicCoords of the inputSequence to determine which SNPs are within
        range. The function searches through all SNPs and creates a bed file of the SNPs which are in range.

        NB: SNP tracks of all SNPs(147) have been downloaded from UCSC for each chromosome
        and stored in the genome/SNPs folder. (GRCh37)
        """
        start = self.genomic_coords[0]
        stop = self.genomic_coords[1]
        if genome_build == 'GRCh37':
            snps = BedTool('/media/sf_S_DRIVE/genomic_resources/primer_design/All-SNPs-GRCh37/' +
                           self.chrom_number + '-GRCh37-All-snps-147.bed.gz')
            # snps = BedTool('/srv/primer_design/s_drive/All-SNPs-GRCh37/' +
            #                self.chrom_number + '-GRCh37-All-snps-147.bed.gz')
        elif genome_build == 'GRCh38':
            snps = BedTool('/media/sf_S_DRIVE/genomic_resources/primer_design/All-SNPs-GRCh38/' +
                           self.chrom_number + '-GRCh38-All-snps-150.bed.gz')
            # snps = BedTool('/srv/primer_design/s_drive/All-SNPs-GRCh38/' +
            #                self.chrom_number + '-GRCh38-All-snps-150.bed.gz')
        target_region = BedTool(
            self.chrom_number + " " + str(start) + " " + str(stop), from_string=True)
        self.snps_bed = snps.intersect(target_region)


class Exon:
    """
    Stores information about the exons. A processed boolean flag is used when determining the target regions.
    Once an exon has been incorporated into a target region, this flag will be set to True.
    """

    def __init__(self, number, loc_start, loc_stop, genom_start, genom_stop):
        self.number = number
        self.loc_start = loc_start
        self.loc_stop = loc_stop
        self.genom_start = genom_start
        self.genom_stop = genom_stop
        self.length = loc_stop - loc_start
        self.processed = False
        self.final_exon = False
        self.primers = []


class TargetRegion:
    """
    Stores information about the target regions. These may be whole exons, multiple exons, or parts of exons depending
    on their size. A masked sequence is created which masks the orginial sequence in the locations where common SNPs
    are found. These will then be avoided by the Primer3 API when designing the primers.
    """

    def __init__(self, target_id, sequence, seq_start, seq_stop, overhang, input_seq):
        self.target_id = target_id
        self.sequence = sequence
        self.seq_start = seq_start
        self.seq_stop = seq_stop
        self.overhang = overhang
        self.input_sequence = input_seq
        self.snps = []
        self.masked_sequence = ""
        self.primers = []

    # def __init__(self, target_id, sequence, seq_start, seq_stop, overhang):
    #     self.target_id = target_id
    #     self.sequence = sequence
    #     self.seq_start = seq_start
    #     self.seq_stop = seq_stop
    #     self.overhang = overhang
    #     self.input_sequence = ""
    #     self.snps = []
    #     self.masked_sequence = ""
    #     self.primers = []

    def set_snps(self):
        """
        Uses the SNP bed file that was created to cover the whole input sequence to search for SNPs which are in
        range of the target region. This saves searching through the whole list of chromosome SNPs for each target
        region. A Snp object is created for each SNP found within range of the target region.
        """
        strand = self.input_sequence.strand
        chromosome_number = self.input_sequence.chrom_number
        input_seq_snps = self.input_sequence.snps_bed
        if strand == "+":
            start = self.seq_start
            stop = self.seq_stop
        else:
            start = self.seq_stop
            stop = self.seq_start
        snps_in_range = []
        target_region = BedTool(chromosome_number + " " + str(start) + " " + str(stop), from_string=True)
        snps_in_region = input_seq_snps.intersect(target_region)
        for snp in snps_in_region:
            try:
                snp_start = int(snp[1])
                snp_stop = int(snp[2])
                snp_id = snp[3]
                ref_ncbi = snp[4]
                observed = snp[5]
                snp_class = snp[6]
                av_het = float(snp[7])
                av_het_se = float(snp[8])
                loc_type = snp[9]
                snp_object = Snp(chromosome_number, snp_start, snp_stop, snp_id, ref_ncbi, observed, snp_class, av_het,
                                 av_het_se, loc_type)
                snps_in_range.append(snp_object)
            except:
                print("error parsing SNP details")
        self.snps = snps_in_range

    def mask_sequence(self, av_het):
        """
        Generates a masked version of the sequence. For each SNP location found within the sequence, the
        nucleotide at that position is changed to an "N"
        """
        # checks to see if the input sequence is on the - strand, and if so converts to the reverse complement
        strand = self.input_sequence.strand
        if strand == "-":
            temp = Seq(self.sequence)
            seq = temp.reverse_complement()
        else:
            seq = self.sequence
        # converts the sequence to a list of characters so the locations of the SNPs can be masked with "N"
        seq_list = list(seq)
        # for each SNP in the region calculates the positional offset from the position of the SNP and
        # the start position of the sequence genomic coordinates
        for snp in self.snps:
            if av_het <= snp.av_het:
                if strand == "+":
                    offset = snp.coord_start + 1 - self.seq_start
                if strand == "-":
                    offset = snp.coord_start - self.seq_start
                # changes the nucleotide at the position of the SNP to an N in the sequence list
                seq_list[offset] = 'N'
        # converts the sequence list back to a string
        seq = ''.join(seq_list)
        # if the sequence was on the - strand the reverse complement is performed again
        if strand == "-":
            temp = Seq(seq)
            seq = temp.reverse_complement()
        # a string interpretation of the masked sequence is assigned to the maskedSequence instance variable
        self.masked_sequence = str(seq)

    def set_primers(self, min_product_size, max_product_size, primer_opt_size, primer_min_size, primer_max_size,
                    primer_opt_tm, primer_min_tm, primer_max_tm, primer_min_gc, primer_max_gc):
        """
        Generates and sets primers for the maskedSequence using Primer3.
        Primers objects are created for each primer and the list of primers is
        assigned to the primers instance variable for inputSequence.

        Currently the Primer3 parameters are fixed.
        """
        # SOP says to ensure primers are positioned at least 30 bases either side of region of interest
        target_start = self.overhang - 30
        target_end = len(self.masked_sequence) - self.overhang + 30
        target_length = target_end - target_start
        results = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': self.target_id,
                'SEQUENCE_TEMPLATE': self.masked_sequence,
                'SEQUENCE_TARGET': [target_start, target_length],
            },
            {
                'PRIMER_PRODUCT_SIZE_RANGE': [min_product_size,max_product_size],
                'PRIMER_OPT_SIZE': primer_opt_size,
                'PRIMER_PICK_LEFT_PRIMER': 1,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_PICK_RIGHT_PRIMER': 1,
                'PRIMER_INTERNAL_MAX_SELF_END': 0,
                'PRIMER_MIN_SIZE': primer_min_size,
                'PRIMER_MAX_SIZE': primer_max_size,
                'PRIMER_INTERNAL_MAX_SIZE': primer_max_size,
                'PRIMER_OPT_TM': primer_opt_tm,
                'PRIMER_MIN_TM': primer_min_tm,
                'PRIMER_MAX_TM': primer_max_tm,
                'PRIMER_MIN_GC': primer_min_gc,
                'PRIMER_MAX_GC': primer_max_gc,
                'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
                'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 0,
                'PRIMER_LIBERAL_BASE': 1,
                'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': 1,
                'PRIMER_LOWERCASE_MASKING': 0,
                'PRIMER_PICK_ANYWAY': 1,
                'PRIMER_EXPLAIN_FLAG': 1,
                'PRIMER_TASK': 'generic',
                'PRIMER_MIN_QUALITY': 0,
                'PRIMER_MIN_END_QUALITY': 0,
                'PRIMER_QUALITY_RANGE_MIN': 0,
                'PRIMER_QUALITY_RANGE_MAX': 100,
                'PRIMER_PAIR_MAX_DIFF_TM': 100,
                'PRIMER_TM_FORMULA': 0,
                'PRIMER_PRODUCT_MIN_TM': -1000000.0,
                'PRIMER_PRODUCT_OPT_TM': 0.0,
                'PRIMER_PRODUCT_MAX_TM': 1000000.0,
                'PRIMER_OPT_GC_PERCENT': 50.0,
                'PRIMER_NUM_RETURN': 5,
                'PRIMER_MAX_END_STABILITY': 9.0,
                'PRIMER_MAX_LIBRARY_MISPRIMING': 12.00,
                'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': 24.00,
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 40.00,
                'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': 70.00,
                'PRIMER_MAX_SELF_ANY_TH': 45.0,
                'PRIMER_MAX_SELF_END_TH': 35.0,
                'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,
                'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
                'PRIMER_MAX_HAIRPIN_TH': 24.0,
                'PRIMER_MAX_TEMPLATE_MISPRIMING': 12.00,
                'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': 24.00,
                'PRIMER_MAX_SELF_ANY': 12.00,
                'PRIMER_MAX_SELF_END': 12.00,
                'PRIMER_PAIR_MAX_COMPL_ANY': 8.00,
                'PRIMER_PAIR_MAX_COMPL_END': 3.00,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_POLY_X': 5,
                'PRIMER_INSIDE_PENALTY': 0,
                'PRIMER_OUTSIDE_PENALTY': 0,
                'PRIMER_GC_CLAMP': 0,
                'PRIMER_MAX_END_GC': 5,
                'PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE': 3,
                'PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE': 3,
                'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION': 7,
                'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 4,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_SALT_CORRECTIONS': 1,
                'PRIMER_SALT_DIVALENT': 0,
                'PRIMER_DNTP_CONC': 0,
                'PRIMER_DNA_CONC': 50.00,
                'PRIMER_SEQUENCING_SPACING': 500,
                'PRIMER_SEQUENCING_INTERVAL': 250,
                'PRIMER_SEQUENCING_LEAD': 50,
                'PRIMER_SEQUENCING_ACCURACY': 20,
                'PRIMER_WT_SIZE_LT': 1.0,
                'PRIMER_WT_SIZE_GT': 1.0,
                'PRIMER_WT_TM_LT': 1.0,
                'PRIMER_WT_TM_GT': 1.0,
                'PRIMER_WT_GC_PERCENT_LT': 0.0,
                'PRIMER_WT_GC_PERCENT_GT': 0.0,
                'PRIMER_WT_SELF_ANY_TH': 0.0,
                'PRIMER_WT_SELF_END_TH': 0.0,
                'PRIMER_WT_HAIRPIN_TH': 0.0,
                'PRIMER_WT_TEMPLATE_MISPRIMING_TH': 0.0,
                'PRIMER_WT_SELF_ANY': 0.0,
                'PRIMER_WT_SELF_END': 0.0,
                'PRIMER_WT_TEMPLATE_MISPRIMING': 0.0,
                'PRIMER_WT_NUM_NS': 0.0,
                'PRIMER_WT_LIBRARY_MISPRIMING': 0.0,
                'PRIMER_WT_SEQ_QUAL': 0.0,
                'PRIMER_WT_END_QUAL': 0.0,
                'PRIMER_WT_POS_PENALTY': 0.0,
                'PRIMER_WT_END_STABILITY': 0.0,
                'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.0,
                'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': 0.0,
                'PRIMER_PAIR_WT_PRODUCT_TM_LT': 0.0,
                'PRIMER_PAIR_WT_PRODUCT_TM_GT': 0.0,
                'PRIMER_PAIR_WT_COMPL_ANY_TH': 0.0,
                'PRIMER_PAIR_WT_COMPL_END_TH': 0.0,
                'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH': 0.0,
                'PRIMER_PAIR_WT_COMPL_ANY': 0.0,
                'PRIMER_PAIR_WT_COMPL_END': 0.0,
                'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING': 0.0,
                'PRIMER_PAIR_WT_DIFF_TM': 0.0,
                'PRIMER_PAIR_WT_LIBRARY_MISPRIMING': 0.0,
                'PRIMER_PAIR_WT_PR_PENALTY': 1.0,
                'PRIMER_PAIR_WT_IO_PENALTY': 0.0,
                'PRIMER_INTERNAL_MIN_SIZE': 18,
                'PRIMER_INTERNAL_OPT_SIZE': 20,
                'PRIMER_INTERNAL_MIN_TM': 57.0,
                'PRIMER_INTERNAL_OPT_TM': 60.0,
                'PRIMER_INTERNAL_MAX_TM': 63.0,
                'PRIMER_INTERNAL_MIN_GC': 20.0,
                'PRIMER_INTERNAL_OPT_GC_PERCENT': 50.0,
                'PRIMER_INTERNAL_MAX_GC': 80.0,
                'PRIMER_INTERNAL_MAX_SELF_ANY_TH': 47.00,
                'PRIMER_INTERNAL_MAX_SELF_END_TH': 47.00,
                'PRIMER_INTERNAL_MAX_HAIRPIN_TH': 47.00,
                'PRIMER_INTERNAL_MAX_SELF_ANY': 12.00,
                'PRIMER_INTERNAL_MIN_QUALITY': 0,
                'PRIMER_INTERNAL_MAX_NS_ACCEPTED': 0,
                'PRIMER_INTERNAL_MAX_POLY_X': 5,
                'PRIMER_INTERNAL_MAX_LIBRARY_MISHYB': 12.00,
                'PRIMER_INTERNAL_SALT_MONOVALENT': 50.0,
                'PRIMER_INTERNAL_DNA_CONC': 50.0,
                'PRIMER_INTERNAL_SALT_DIVALENT': 0.0,
                'PRIMER_INTERNAL_DNTP_CONC': 0.0,
                'PRIMER_INTERNAL_WT_SIZE_LT': 1.0,
                'PRIMER_INTERNAL_WT_SIZE_GT': 1.0,
                'PRIMER_INTERNAL_WT_TM_LT': 1.0,
                'PRIMER_INTERNAL_WT_TM_GT': 1.0,
                'PRIMER_INTERNAL_WT_GC_PERCENT_LT': 0.0,
                'PRIMER_INTERNAL_WT_GC_PERCENT_GT': 0.0,
                'PRIMER_INTERNAL_WT_SELF_ANY_TH': 0.0,
                'PRIMER_INTERNAL_WT_SELF_END_TH': 0.0,
                'PRIMER_INTERNAL_WT_HAIRPIN_TH': 0.0,
                'PRIMER_INTERNAL_WT_SELF_ANY': 0.0,
                'PRIMER_INTERNAL_WT_SELF_END': 0.0,
                'PRIMER_INTERNAL_WT_NUM_NS': 0.0,
                'PRIMER_INTERNAL_WT_LIBRARY_MISHYB': 0.0,
                'PRIMER_INTERNAL_WT_SEQ_QUAL': 0.0,
                'PRIMER_INTERNAL_WT_END_QUAL': 0.0,
            })
        num_pair_primers = results.get('PRIMER_PAIR_NUM_RETURNED')
        # creates a list containing the required number of primer objects
        primers = [Primer() for _ in range(num_pair_primers)]
        count = 0
        # The output of Primer3 is parsed to extract the primer information which is assigned to the
        # Primer object instance variables
        for primer in primers:
            primer.product_size = results.get('PRIMER_PAIR_' + str(count) + '_PRODUCT_SIZE')
            primer.forward_seq = results.get('PRIMER_LEFT_' + str(count) + '_SEQUENCE')
            primer.forward_start = results.get('PRIMER_LEFT_' + str(count))[0]
            primer.forward_length = results.get('PRIMER_LEFT_' + str(count))[1]
            primer.forward_tm = results.get('PRIMER_LEFT_' + str(count) + '_TM')
            primer.forward_gc = results.get('PRIMER_LEFT_' + str(count) + '_GC_PERCENT')
            primer.reverse_seq = results.get('PRIMER_RIGHT_' + str(count) + '_SEQUENCE')
            primer.reverse_start = results.get('PRIMER_RIGHT_' + str(count))[0]
            primer.reverse_length = results.get('PRIMER_RIGHT_' + str(count))[1]
            primer.reverse_tm = results.get('PRIMER_RIGHT_' + str(count) + '_TM')
            primer.reverse_gc = results.get('PRIMER_RIGHT_' + str(count) + '_GC_PERCENT')
            primer.internal_seq = results.get('PRIMER_INTERNAL_' + str(count) + '_SEQUENCE')
            primer.internal_start = results.get('PRIMER_INTERNAL_' + str(count))[0]
            primer.internal_length = results.get('PRIMER_INTERNAL_' + str(count))[1]
            primer.internal_tm = results.get('PRIMER_INTERNAL_' + str(count) + '_TM')
            primer.internal_gc = results.get('PRIMER_INTERNAL_' + str(count) + '_GC_PERCENT')
            if self.input_sequence.strand == "+":
                primer.forward_genomic_coords = (self.seq_start + primer.forward_start,
                                                 self.seq_start + primer.forward_start + primer.forward_length - 1)
                primer.reverse_genomic_coords = (self.seq_start + primer.reverse_start - primer.reverse_length + 1,
                                                 self.seq_start + primer.reverse_start)
                primer.internal_genomic_coords = (self.seq_start + primer.internal_start,
                                                  self.seq_start + primer.internal_start + primer.internal_length - 1)
            elif self.input_sequence.strand == "-":
                primer.forward_genomic_coords = (self.seq_start - primer.forward_start - primer.forward_length + 1,
                                                 self.seq_start - primer.forward_start)
                primer.reverse_genomic_coords = (self.seq_start - primer.reverse_start,
                                                 self.seq_start - primer.reverse_start + primer.reverse_length - 1)
                primer.internal_genomic_coords = (self.seq_start - primer.internal_start,
                                                  self.seq_start - primer.internal_start - primer.internal_length + 1)
            count += 1
        self.primers = primers

    def set_bespoke_primers(self, min_product_size, max_product_size, primer_opt_size, primer_min_size, primer_max_size,
                    primer_opt_tm, primer_min_tm, primer_max_tm, primer_min_gc, primer_max_gc):
        """
        Generates and sets primers for the maskedSequence using Primer3.
        Primers objects are created for each primer and the list of primers is
        assigned to the primers instance variable for inputSequence.

        Currently the Primer3 parameters are fixed.
        """
        target_start = self.overhang #- 5
        target_end = len(self.masked_sequence) - self.overhang #+ 5
        target_length = target_end - target_start
        results = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': self.target_id,
                'SEQUENCE_TEMPLATE': self.masked_sequence,
                'SEQUENCE_TARGET': [target_start, target_length],
            },
            {
                'PRIMER_PRODUCT_SIZE_RANGE': [min_product_size,max_product_size],
                'PRIMER_OPT_SIZE': primer_opt_size,
                'PRIMER_PICK_LEFT_PRIMER': 1,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_PICK_RIGHT_PRIMER': 1,
                'PRIMER_INTERNAL_MAX_SELF_END': 0,
                'PRIMER_MIN_SIZE': primer_min_size,
                'PRIMER_MAX_SIZE': primer_max_size,
                'PRIMER_INTERNAL_MAX_SIZE': primer_max_size,
                'PRIMER_OPT_TM': primer_opt_tm,
                'PRIMER_MIN_TM': primer_min_tm,
                'PRIMER_MAX_TM': primer_max_tm,
                'PRIMER_MIN_GC': primer_min_gc,
                'PRIMER_MAX_GC': primer_max_gc,
                'PRIMER_LIBERAL_BASE': 1,
                'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': 1,
                'PRIMER_LOWERCASE_MASKING': 0,
                'PRIMER_MIN_QUALITY': 0,
                'PRIMER_MIN_END_QUALITY': 0,
                'PRIMER_QUALITY_RANGE_MIN': 0,
                'PRIMER_QUALITY_RANGE_MAX': 100,
                'PRIMER_PAIR_MAX_DIFF_TM': 100,
                'PRIMER_TM_FORMULA': 0,
                'PRIMER_NUM_RETURN': 5,
                'PRIMER_MAX_END_STABILITY': 9.0,
                'PRIMER_MAX_LIBRARY_MISPRIMING': 12.00,
                'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': 24.00,
                'PRIMER_MAX_TEMPLATE_MISPRIMING': 12.00,
                'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': 24.00,
                'PRIMER_MAX_SELF_ANY': 12.00,
                'PRIMER_MAX_SELF_END': 12.00,
                'PRIMER_PAIR_MAX_COMPL_ANY': 8.00,
                'PRIMER_PAIR_MAX_COMPL_END': 3.00,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_POLY_X': 5,
                'PRIMER_INSIDE_PENALTY': 0,
                'PRIMER_OUTSIDE_PENALTY': 0,
                'PRIMER_GC_CLAMP': 0,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_SALT_CORRECTIONS': 1,
                'PRIMER_SALT_DIVALENT': 0,
                'PRIMER_DNTP_CONC': 0,
                'PRIMER_DNA_CONC': 50.00,
                'PRIMER_WT_SIZE_LT': 1.0,
                'PRIMER_WT_SIZE_GT': 1.0,
                'PRIMER_WT_TM_LT': 1.0,
                'PRIMER_WT_TM_GT': 1.0,
                'PRIMER_WT_GC_PERCENT_LT': 0.0,
                'PRIMER_WT_GC_PERCENT_GT': 0.0,
                'PRIMER_WT_SELF_ANY': 0.0,
                'PRIMER_WT_SELF_END': 0.0,
                'PRIMER_WT_TEMPLATE_MISPRIMING': 0.0,
                'PRIMER_WT_NUM_NS': 0.0,
                'PRIMER_WT_LIBRARY_MISPRIMING': 0.0,
                'PRIMER_WT_SEQ_QUAL': 0.0,
                'PRIMER_WT_END_QUAL': 0.0,
                'PRIMER_WT_POS_PENALTY': 0.0,
                'PRIMER_WT_END_STABILITY': 0.0,
                'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.0,
                'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': 0.0,
                'PRIMER_PAIR_WT_PRODUCT_TM_LT': 0.0,
                'PRIMER_PAIR_WT_PRODUCT_TM_GT': 0.0,
                'PRIMER_PAIR_WT_COMPL_ANY': 0.0,
                'PRIMER_PAIR_WT_COMPL_END': 0.0,
                'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING': 0.0,
                'PRIMER_PAIR_WT_DIFF_TM': 0.0,
                'PRIMER_PAIR_WT_LIBRARY_MISPRIMING': 0.0,
                'PRIMER_PAIR_WT_PR_PENALTY': 1.0,
                'PRIMER_PAIR_WT_IO_PENALTY': 0.0,
                'PRIMER_INTERNAL_MIN_SIZE': 18,
                'PRIMER_INTERNAL_OPT_SIZE': 20,
                'PRIMER_INTERNAL_MIN_TM': 57.0,
                'PRIMER_INTERNAL_OPT_TM': 60.0,
                'PRIMER_INTERNAL_MAX_TM': 63.0,
                'PRIMER_INTERNAL_MIN_GC': 20.0,
                'PRIMER_INTERNAL_MAX_GC': 80.0,
                'PRIMER_INTERNAL_MAX_SELF_ANY': 12.00,
                'PRIMER_INTERNAL_MIN_QUALITY': 0,
                'PRIMER_INTERNAL_MAX_NS_ACCEPTED': 0,
                'PRIMER_INTERNAL_MAX_POLY_X': 5,
                'PRIMER_INTERNAL_MAX_LIBRARY_MISHYB': 12.00,
                'PRIMER_INTERNAL_SALT_MONOVALENT': 50.0,
                'PRIMER_INTERNAL_DNA_CONC': 50.0,
                'PRIMER_INTERNAL_SALT_DIVALENT': 0.0,
                'PRIMER_INTERNAL_DNTP_CONC': 0.0,
                'PRIMER_INTERNAL_WT_SIZE_LT': 1.0,
                'PRIMER_INTERNAL_WT_SIZE_GT': 1.0,
                'PRIMER_INTERNAL_WT_TM_LT': 1.0,
                'PRIMER_INTERNAL_WT_TM_GT': 1.0,
                'PRIMER_INTERNAL_WT_GC_PERCENT_LT': 0.0,
                'PRIMER_INTERNAL_WT_GC_PERCENT_GT': 0.0,
                'PRIMER_INTERNAL_WT_SELF_ANY': 0.0,
                'PRIMER_INTERNAL_WT_NUM_NS': 0.0,
                'PRIMER_INTERNAL_WT_LIBRARY_MISHYB': 0.0,
                'PRIMER_INTERNAL_WT_SEQ_QUAL': 0.0,
            })
        num_pair_primers = results.get('PRIMER_PAIR_NUM_RETURNED')
        # creates a list containing the required number of primer objects
        primers = [Primer() for _ in range(num_pair_primers)]
        count = 0
        # The output of Primer3 is parsed to extract the primer information which is assigned to the
        # Primer object instance variables
        for primer in primers:
            primer.product_size = results.get('PRIMER_PAIR_' + str(count) + '_PRODUCT_SIZE')
            primer.forward_seq = results.get('PRIMER_LEFT_' + str(count) + '_SEQUENCE')
            primer.forward_start = results.get('PRIMER_LEFT_' + str(count))[0]
            primer.forward_length = results.get('PRIMER_LEFT_' + str(count))[1]
            primer.forward_tm = results.get('PRIMER_LEFT_' + str(count) + '_TM')
            primer.forward_gc = results.get('PRIMER_LEFT_' + str(count) + '_GC_PERCENT')
            primer.reverse_seq = results.get('PRIMER_RIGHT_' + str(count) + '_SEQUENCE')
            primer.reverse_start = results.get('PRIMER_RIGHT_' + str(count))[0]
            primer.reverse_length = results.get('PRIMER_RIGHT_' + str(count))[1]
            primer.reverse_tm = results.get('PRIMER_RIGHT_' + str(count) + '_TM')
            primer.reverse_gc = results.get('PRIMER_RIGHT_' + str(count) + '_GC_PERCENT')
            
            if self.input_sequence.strand == "+":
                primer.forward_genomic_coords = (self.seq_start + primer.forward_start,
                                                 self.seq_start + primer.forward_start + primer.forward_length - 1)
                primer.reverse_genomic_coords = (self.seq_start + primer.reverse_start - primer.reverse_length + 1,
                                                 self.seq_start + primer.reverse_start)
                primer.internal_genomic_coords = (self.seq_start + primer.internal_start,
                                                  self.seq_start + primer.internal_start + primer.internal_length - 1)
            elif self.input_sequence.strand == "-":
                primer.forward_genomic_coords = (self.seq_start - primer.forward_start - primer.forward_length + 1,
                                                 self.seq_start - primer.forward_start)
                primer.reverse_genomic_coords = (self.seq_start - primer.reverse_start,
                                                 self.seq_start - primer.reverse_start + primer.reverse_length - 1)
                primer.internal_genomic_coords = (self.seq_start - primer.internal_start,
                                                  self.seq_start - primer.internal_start - primer.internal_length + 1)
            count += 1
        self.primers = primers

class Primer:
    """
    Primer class for holding information relating to a set of primers (i.e. forward, reverse and internal)
    """

    def __init__(self):
        self.product_size = 0
        self.forward_seq = ""
        self.forward_start = 0
        self.forward_length = 0
        self.forward_tm = 0.0
        self.forward_gc = 0.0
        self.forward_genomic_coords = None
        self.forward_snps = []
        self.reverse_seq = ""
        self.reverse_start = 0
        self.reverse_length = 0
        self.reverse_tm = 0.0
        self.reverse_gc = 0.0
        self.reverse_genomic_coords = None
        self.reverse_snps = []
        self.internal_seq = ""
        self.internal_start = 0
        self.internal_length = 0
        self.internal_tm = 0.0
        self.internal_gc = 0.0
        self.internal_genomic_coords = None
        self.internal_snps = []

    def set_snps(self, target_region):
        """
        Sets the SNPs that are found within the region of the primers
        """
        self.forward_snps = self.find_snps(self.forward_genomic_coords, target_region)
        self.reverse_snps = self.find_snps(self.reverse_genomic_coords, target_region)
        self.internal_snps = self.find_snps(self.internal_genomic_coords, target_region)

    @staticmethod
    def find_snps(genomic_coords, target_region):
        """
        Helper method to find SNPs location within the region of the primer
        """
        target_region_snps = target_region.snps
        primer_site_snps = []
        for snp in target_region_snps:
            if (genomic_coords[0] - 1 <= snp.coord_start < genomic_coords[1]) or (
                            genomic_coords[0] - 1  < snp.coord_end <= genomic_coords[1]):
                primer_site_snps.append(snp)
        return primer_site_snps


class Snp:
    """
    Snp class holds information about individial SNPs.
    """

    def __init__(self, chrom_number, snp_start, snp_stop, snp_id, ref_ncbi, observed, snp_class, av_het, av_het_se,
                 loc_type):
        self.chrom_number = chrom_number
        self.coord_start = snp_start
        self.coord_end = snp_stop
        self.snp_id = snp_id
        self.ref_NCBI = ref_ncbi
        self.observed = observed
        self.snp_class = snp_class
        self.av_het = av_het
        self.av_het_se = av_het_se
        self.loc_type = loc_type


# main functions. 
def get_ref_seq(ng_id):
    """
    Gets the reference sequence from NCBI using the NG_ number supplied by the user.
    """
    Entrez.email = "stone_edward@hotmail.com"
    handle = Entrez.efetch(db="nucleotide", id=ng_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    return record


def write_to_csv(input_sequence):
    """
    Writes primers that have been designed to a CSV file
    """
    # inputSequence ID used as the filename
    filename = input_sequence.gene_name + '-' + input_sequence.refseq_id
    targets = input_sequence.target_regions
    if input_sequence.strand == '+':
        strand = "Forward (+)"
    else:
        strand = "Reverse (-)"
    # with open('output/' + filename + '.csv', 'w') as csvfile:

    # with open('/srv/primer_design/s_drive/designs/' + input_sequence.gene_name + '.csv', 'w') as csvfile:
    with open('/media/sf_S_DRIVE/genomic_resources/primer_design/designs/' + input_sequence.gene_name + '.csv', 'w') as csvfile:
        f = csv.writer(csvfile, delimiter=',',
                       quotechar=',', quoting=csv.QUOTE_MINIMAL)
        # f.writerow(['Gene:', filename, '\n'])
        # f.writerow(['Designed primers:'])
        f.writerow(['Gene', 'NG number', 'Strand', 'Target (Exon)', 'Product size', 'Forward primer sequence', 'Genomic Coords', 'Forward TM',
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
                f.writerow([input_sequence.gene_name, input_sequence.refseq_id, strand, target.target_id, primer.product_size, primer.forward_seq,
                            input_sequence.chrom_number + ":" + str(primer.forward_genomic_coords[0]) + "-" + str(
                                primer.forward_genomic_coords[1]), round(primer.forward_tm, 2),
                            round(primer.forward_gc, 2), forward_snps, primer.reverse_seq,
                            input_sequence.chrom_number + ":" + str(primer.reverse_genomic_coords[0]) + "-" + str(
                                primer.reverse_genomic_coords[1]),
                            round(primer.reverse_tm, 2), round(primer.reverse_gc, 2), reverse_snps])


def write_to_bed(input_sequence):
    """
    Writes primers that have been designed to a CSV file
    """
    # inputSequence ID used as the filename
    filename = input_sequence.gene_name + '-' + input_sequence.refseq_id
    targets = input_sequence.target_regions
    # with open('/srv/primer_design/s_drive/designs/' + input_sequence.gene_name + '.bed', 'w') as csvfile:
    with open('/media/sf_S_DRIVE/genomic_resources/primer_design/designs/' + input_sequence.gene_name + '.bed', 'w') as csvfile:
        f = csv.writer(csvfile, delimiter='\t',
                       quotechar=';', quoting=csv.QUOTE_MINIMAL)
        f.writerow(['track name="' + filename + '" description=' + '"Primers designed for' + filename +
                    '" visibility=2 itemRgb="On"'])
        if input_sequence.strand == '+':
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

        else:
            for target in targets:
                f.writerow([input_sequence.chrom_number, target.seq_stop + target.overhang - 30,
                            target.seq_start - target.overhang + 30, target.target_id, 0, input_sequence.strand,
                            target.seq_stop + target.overhang - 30, target.seq_start - target.overhang + 30,
                            '255,0,0'])
            for target in targets:
                for primer in target.primers:
                    f.writerow([input_sequence.chrom_number, primer.reverse_genomic_coords[0] - 1,
                                primer.forward_genomic_coords[1], target.target_id, 0, input_sequence.strand,
                                primer.reverse_genomic_coords[1], primer.forward_genomic_coords[0] - 1, '0,0,255'])


def design_by_symbol(genome_build, ng_number, chromosome, start, strand, max_avhet, min_prod_size,
        max_product_size, primer_opt_size, primer_min_size, primer_max_size, primer_opt_tm, primer_min_tm,
        primer_max_tm, primer_min_gc, primer_max_gc):
    """
    Main function for script.
    Accepts a ref_seq NG_ id as input and calls functions to generate primers for this input.
    """
    ref_seq = get_ref_seq(ng_number)
    # max produce size minus 30b either side of region and at least 100b either side for primer binding as stated
    # in SOP (total = -260)
    target_region_max_size = max_product_size - 260
    if int(min_prod_size / 2) > 250:
        overhang = int(min_prod_size / 2)
    else:
        # size of overhang sequence available for primers to bind to around target sites set at least 250 bases
        overhang = 250
    # Creates an input_sequence object using the record ID and sequence from FASTA input.
    input_sequence = InputSequence(ref_seq.id, str(ref_seq.seq), chrom_number=chromosome,
                                   genomic_coords=(start,start+len(str(ref_seq.seq))-1), strand=strand)
    # Performs a BLAT search of the input sequenece to determine genomic location.
    #input_sequence.blat_search(ref_seq.format("fasta"), genome_build)
    # Setter methods called to calculate and assign other instance variables.
    #input_sequence.set_genomic_location(psl_file_path='output/gfOutput.psl')
    input_sequence.set_exons(ref_seq)
    input_sequence.set_snps_bed(genome_build)
    input_sequence.set_target_regions(target_region_max_size, overhang)
    target_regions = input_sequence.target_regions
    for target in target_regions:
        target.set_snps()
        target.mask_sequence(max_avhet)
        target.set_primers(min_prod_size, max_product_size, primer_opt_size, primer_min_size, primer_max_size,
                           primer_opt_tm, primer_min_tm, primer_max_tm, primer_min_gc, primer_max_gc)
        primers = target.primers
        for primer in primers:
            primer.set_snps(target)
    # Primers that have been designed are output to a CSV file.
    write_to_csv(input_sequence)
    write_to_bed(input_sequence)
    return input_sequence

if __name__ == '__main__':
    gene_symbol = sys.argv[1]
    print('gene symbol = ', gene_symbol)
    type(gene_symbol)
    main(gene_symbol, 1, 300, 750, 20, 18, 27, 60, 57, 63, 20, 80)


