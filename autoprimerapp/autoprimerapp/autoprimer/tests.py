import os
from django.test import TestCase
from autoprimer import autoprimer
from urllib.error import HTTPError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pybedtools import BedTool


class AutoprimerTestCase(TestCase):
    def setUp(self):
        # dummy InputSequence object with sequence length = 400
        self.test_input_seq = autoprimer.InputSequence('NG_TEST01',
                                                       'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                       'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                       'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                       'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                       'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                       'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                       'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                       'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC')

        # set up dummy RefSeq object
        f1 = SeqFeature(FeatureLocation(5, 10, strand=1), type="exon", qualifiers={'gene': "TESTG", 'number': [1]})
        f2 = SeqFeature(FeatureLocation(13, 25, strand=1), type="exon", qualifiers={'gene': "TESTG", 'number': [2]})
        f3 = SeqFeature(FeatureLocation(31, 40, strand=1), type="exon", qualifiers={'gene': "TESTG", 'number': [3]})
        feature_list = [f1, f2, f3]

        self.test_ref_seq = SeqRecord(Seq("AGTCATCTACATCTACTGTGCATCAGCTAGCTGTCGATCGATCGATCGTACGTCGATCGT",
                                          IUPAC.unambiguous_dna),
                                      id="testRefSeq1", name="Test", description="a test file", features=feature_list)

        # dummy Exon objects created
        test_exon_1 = autoprimer.Exon(1, 10, 20, 10, 20)
        test_exon_2 = autoprimer.Exon(2, 30, 40, 30, 40)
        test_exon_3 = autoprimer.Exon(3, 100, 200, 100, 200)
        test_exon_4 = autoprimer.Exon(4, 350, 360, 350, 360)
        test_exon_4.final_exon = True
        self.test_exons_positive_strand = [test_exon_1, test_exon_2, test_exon_3, test_exon_4]
        test_exon_5 = autoprimer.Exon(1, 10, 20, 390, 380)
        test_exon_6 = autoprimer.Exon(2, 30, 40, 370, 360)
        test_exon_7 = autoprimer.Exon(3, 100, 200, 300, 200)
        test_exon_8 = autoprimer.Exon(4, 350, 360, 50, 40)
        test_exon_8.final_exon = True
        self.test_exons_negative_strand = [test_exon_5, test_exon_6, test_exon_7, test_exon_8]

        # dummy TargetRegion objects created
        self.test_target_region = autoprimer.TargetRegion(1, 'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                          'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                          'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                                                          'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC',
                                                          0, 200, 60, self.test_input_seq)

        # dummy SNP bed files
        snp_string = """
        chr1 0 1 rsTESTsnp1 A A/T single 0.02 0.001 exact
        chr1 49 50 rsTESTsnp2 C C/G single 0.004 0.001 exact
        chr1 199 200 rsTESTsnp3 C C/A single 0.06 0.001 exact
        chr1 200 201 rsTESTsnp4 A A/G single 0.003 0.001 exact
        chr1 248 249 rsTESTsnp5 T T/C single 0.3 0.001 exact
        """
        self.test_input_seq_snp_file = BedTool(snp_string, from_string=True)

        # dummy SNP objects created
        snp_1 = autoprimer.Snp('chr1', 0, 1, "rsTESTsnp1", "A", "A/T", "single", 0.02, 0.001, "exact")
        snp_2 = autoprimer.Snp('chr1', 49, 50, "rsTESTsnp2", "C", "C/G", "single", 0.004, 0.001, "exact")
        snp_3 = autoprimer.Snp('chr1', 199, 200, "rsTESTsnp3", "C", "C/A", "single", 0.06, 0.001, "exact")
        self.test_snp_list = [snp_1, snp_2, snp_3]
        snp_4 = autoprimer.Snp('chr1', 9, 10, "rsTESTsnp4", "C", "C/T", "single", 0.002, 0.001, "exact")
        snp_5 = autoprimer.Snp('chr1', 10, 11, "rsTESTsnp5", "A", "A/G", "single", 0.004, 0.001, "exact")
        snp_6 = autoprimer.Snp('chr1', 15, 16, "rsTESTsnp6", "T", "T/A", "single", 0.006, 0.001, "exact")
        snp_7 = autoprimer.Snp('chr1', 29, 30, "rsTESTsnp7", "C", "C/T", "single", 0.002, 0.001, "exact")
        snp_8 = autoprimer.Snp('chr1', 30, 31, "rsTESTsnp8", "A", "A/G", "single", 0.004, 0.001, "exact")
        snp_9 = autoprimer.Snp('chr1', 169, 170, "rsTESTsnp9", "A", "A/T", "single", 0.02, 0.001, "exact")
        snp_10 = autoprimer.Snp('chr1', 170, 171, "rsTESTsnp10", "G", "G/A", "single", 0.004, 0.001, "exact")
        snp_11 = autoprimer.Snp('chr1', 180, 181, "rsTESTsnp11", "G", "G/A", "single", 0.006, 0.001, "exact")
        snp_12 = autoprimer.Snp('chr1', 189, 190, "rsTESTsnp12", "A", "A/T", "single", 0.002, 0.001, "exact")
        snp_13 = autoprimer.Snp('chr1', 190, 191, "rsTESTsnp13", "G", "G/T", "single", 0.004, 0.001, "exact")
        self.test_primer_snp_list = [snp_4, snp_5, snp_6, snp_7, snp_8, snp_9, snp_10, snp_11, snp_12, snp_13]

        # dummy Primer object created
        self.test_primer = autoprimer.Primer()
        self.test_primer.forward_start = 10
        self.test_primer.forward_length = 20
        self.test_primer.forward_genomic_coords = (10, 30)
        self.test_primer.reverse_start = 190
        self.test_primer.reverse_length = 20
        self.test_primer.reverse_genomic_coords = (170, 190)
        self.test_primer.internal_genomic_coords = (100, 110)

    def test_get_ref_seq_return_type(self):
        ref_seq = autoprimer.get_ref_seq('NG_009896')
        # checks that the get_ref_seq() function is returning the expected object type
        self.assertEqual(str(type(ref_seq)), "<class 'Bio.SeqRecord.SeqRecord'>")

    def test_get_ref_seq_correct_record_returned(self):
        ref_seq = autoprimer.get_ref_seq('NG_009896')
        # checks that the ID of the returned object is correct
        self.assertIn('NG_009896', ref_seq.id)

    def test_get_ref_seq_error_thrown_with_invalid_entry(self):
        # checks that an HTTP error is raised when input is invalid
        with self.assertRaises(HTTPError):
            autoprimer.get_ref_seq('Bob')

    # def test_blat_search_psl_file_created(self):
    #     # checks that a results file is created following the blat search
    #     if os.path.exists('output/gfOutput.psl'):
    #         os.remove('output/gfOutput.psl')
    #     self.assertFalse(os.path.exists('output/gfOutput.psl'))
    #     f = open('testdata/fasta.txt', 'r')
    #     autoprimer.InputSequence.blat_search(f.read())
    #     self.assertTrue(os.path.exists('output/gfOutput.psl'))

    def test_set_genomic_location_successful_hit_found(self):
        # checks that the genomic location variables are set correctly following a successful match
        self.test_input_seq.sequence_length = 29803
        self.test_input_seq.set_genomic_location('testdata/gfOutputTestSuccess.psl')
        self.assertEqual(self.test_input_seq.chrom_number, 'chrX')
        self.assertEqual(self.test_input_seq.genomic_coords, (153765459, 153795261))
        self.assertEqual(self.test_input_seq.strand, '+')

    def test_set_genomic_location_no_successful_hit_found(self):
        # checks that a ValueError is raised when no match is found
        with self.assertRaises(ValueError):
            self.test_input_seq.set_genomic_location('testdata/gfOutputTestFail.psl')

    def test_set_exons_positive_strand(self):
        self.test_input_seq.strand = "+"
        self.test_input_seq.genomic_coords = (100, 200)
        self.test_input_seq.set_exons(self.test_ref_seq)
        count = 0
        for exon in self.test_input_seq.exons:
            if exon.number == 1:
                self.assertEqual(exon.loc_start, 4)
                self.assertEqual(exon.loc_stop, 9)
                self.assertEqual(exon.genom_start, 104)
                self.assertEqual(exon.genom_stop, 109)
                self.assertFalse(exon.final_exon)
                count += 1
            if exon.number == 2:
                self.assertEqual(exon.loc_start, 12)
                self.assertEqual(exon.loc_stop, 24)
                self.assertEqual(exon.genom_start, 112)
                self.assertEqual(exon.genom_stop, 124)
                self.assertFalse(exon.final_exon)
                count += 1
            if exon.number == 3:
                self.assertEqual(exon.loc_start, 30)
                self.assertEqual(exon.loc_stop, 39)
                self.assertEqual(exon.genom_start, 130)
                self.assertEqual(exon.genom_stop, 139)
                self.assertTrue(exon.final_exon)
                count += 1
        # checks that all expected exons exist and have been checked
        self.assertEqual(count, 3)

    def test_set_exons_negative_strand(self):
        self.test_input_seq.strand = "-"
        self.test_input_seq.genomic_coords = (100, 200)
        self.test_input_seq.set_exons(self.test_ref_seq)
        count = 0
        for exon in self.test_input_seq.exons:
            if exon.number == 1:
                self.assertEqual(exon.loc_start, 6)
                self.assertEqual(exon.loc_stop, 11)
                self.assertEqual(exon.genom_start, 194)
                self.assertEqual(exon.genom_stop, 189)
                self.assertFalse(exon.final_exon)
                count += 1
            if exon.number == 2:
                self.assertEqual(exon.loc_start, 14)
                self.assertEqual(exon.loc_stop, 26)
                self.assertEqual(exon.genom_start, 186)
                self.assertEqual(exon.genom_stop, 174)
                self.assertFalse(exon.final_exon)
                count += 1
            if exon.number == 3:
                self.assertEqual(exon.loc_start, 32)
                self.assertEqual(exon.loc_stop, 41)
                self.assertEqual(exon.genom_start, 168)
                self.assertEqual(exon.genom_stop, 159)
                self.assertTrue(exon.final_exon)
                count += 1
        # checks that all expected exons exist and have been checked
        self.assertEqual(count, 3)

    def test_set_snps_bed(self):
        # All SNP version 147 has 3 SNPs from chr22:16050035-16050114. This includes 1 SNP at each extreme.
        self.test_input_seq.chrom_number = "chr22"
        self.test_input_seq.genomic_coords = (16050035, 16050115)
        self.test_input_seq.set_snps_bed('GRCh37')
        # Check that 3 SNPs are recorded within range
        self.assertEqual(len(self.test_input_seq.snps_bed), 3)

    def test_set_target_regions_positive_strand(self):
        self.test_input_seq.exons = self.test_exons_positive_strand
        self.test_input_seq.strand = "+"
        self.test_input_seq.set_target_regions(30, 10)
        count = 0
        for target in self.test_input_seq.target_regions:
            if target.target_id == "1-2":
                self.assertEqual(target.seq_start, 0)
                self.assertEqual(target.seq_stop, 50)
                count += 1
            if target.target_id == "3.1":
                self.assertEqual(target.seq_start, 90)
                self.assertEqual(target.seq_stop, 135)
                count += 1
            if target.target_id == "3.2":
                self.assertEqual(target.seq_start, 115)
                self.assertEqual(target.seq_stop, 160)
                count += 1
            if target.target_id == "3.3":
                self.assertEqual(target.seq_start, 140)
                self.assertEqual(target.seq_stop, 185)
                count += 1
            if target.target_id == "3.4":
                self.assertEqual(target.seq_start, 165)
                self.assertEqual(target.seq_stop, 210)
                count += 1
            if target.target_id == 4:
                self.assertEqual(target.seq_start, 340)
                self.assertEqual(target.seq_stop, 370)
                count += 1
        # checks that each expected target region exists and has been checked
        self.assertEqual(len(self.test_input_seq.target_regions), count)

    def test_set_target_regions_negative_strand(self):
        self.test_input_seq.exons = self.test_exons_negative_strand
        self.test_input_seq.strand = "-"
        self.test_input_seq.set_target_regions(30, 10)
        count = 0
        for target in self.test_input_seq.target_regions:
            if target.target_id == "1-2":
                self.assertEqual(target.seq_stop, 350)
                self.assertEqual(target.seq_start, 400)
                count += 1
            if target.target_id == "3.1":
                self.assertEqual(target.seq_stop, 265)
                self.assertEqual(target.seq_start, 310)
                count += 1
            if target.target_id == "3.2":
                self.assertEqual(target.seq_stop, 240)
                self.assertEqual(target.seq_start, 285)
                count += 1
            if target.target_id == "3.3":
                self.assertEqual(target.seq_stop, 215)
                self.assertEqual(target.seq_start, 260)
                count += 1
            if target.target_id == "3.4":
                self.assertEqual(target.seq_stop, 190)
                self.assertEqual(target.seq_start, 235)
                count += 1
            if target.target_id == 4:
                self.assertEqual(target.seq_stop, 30)
                self.assertEqual(target.seq_start, 60)
                count += 1
        # checks that each expected target region exists and has been checked
        self.assertEqual(len(self.test_input_seq.target_regions), count)

    def test_target_set_snps(self):
        self.test_input_seq.snps_bed = self.test_input_seq_snp_file
        self.test_input_seq.chrom_number = 'chr1'
        self.test_input_seq.strand = "+"
        self.test_target_region.set_snps()
        self.assertEqual(len(self.test_target_region.snps), 3)

    def test_mask_sequence_positve_strand_mask_all(self):
        self.test_input_seq.strand = "+"
        self.test_target_region.snps = self.test_snp_list
        self.test_target_region.mask_sequence(0.00)
        self.assertEqual(self.test_target_region.masked_sequence, 'NGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTN' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTN')

    def test_mask_sequence_positve_strand_mask_over_0_02(self):
        self.test_input_seq.strand = "+"
        self.test_target_region.snps = self.test_snp_list
        self.test_target_region.mask_sequence(0.02)
        self.assertEqual(self.test_target_region.masked_sequence, 'NGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTN')

    def test_mask_sequence_negative_strand_mask_all(self):
        self.test_input_seq.strand = "-"
        self.test_target_region.snps = self.test_snp_list
        self.test_target_region.mask_sequence(0.00)
        self.assertEqual(self.test_target_region.masked_sequence, 'NGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'NGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTN')

    def test_mask_sequence_negative_strand_mask_over_0_03(self):
        self.test_input_seq.strand = "-"
        self.test_target_region.snps = self.test_snp_list
        self.test_target_region.mask_sequence(0.03)
        self.assertEqual(self.test_target_region.masked_sequence, 'NGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC' +
                         'AGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTCAGTCATGGTC')

    def test_set_primers(self):
        self.test_target_region.masked_sequence = ('GGCGCCCCCGAGCCGCCCTGCTGCGCCGCGCCCGACGCCGCAGCCGCAGCCTTCCCGCCCTGCGC' +
                                                   'TGCCGCCGCCTCCCCGCCACTCTACTCGCAGGTCCCCGACCGCCTGGTACTGCCCGCGACGCGCC' +
                                                   'CCGGCCCCGGCCCGCTGCCCGCTGAGCCCCNNCTGGCCTTGGCCGGGCCGGCAGCCGCTCTCGGC' +
                                                   'CCGCTCAGCCCTGGGGAGGCCTACCTGAGGCAGCCGGGCTTCGCGTCGGGGCTGGAGCGCTACCT' +
                                                   'GTGAGCCTGCGCCGCGCGGGCAGGCACCTGTGCGACCTGTGCCCCGGACCTGCGGCGCCGCCCTC' +
                                                   'GAGCGCCCCNTCTCTACCCCCCACCCTGGCTTGGAGCACACCCTGCGCCTCTCGTCGTGGCCTCC' +
                                                   'TGGACTCAAACTCCTGCTACATCCTTCTCCGTCCCCCATCCCTGGGGAGGCTCCCACCATCTCGC' +
                                                   'CACTGGACAGAAGCGTCCCCTTTGACCTGCCAGCCTCTCATTTCTTCTCCCTCCACTTGTGAGCG' +
                                                   'CCCCCAGCTTTGCGGCGCCCCCCACCCCAGCGCTGTCTGTGGGTCCCTTGCCCCGGAGCTGACTC' +
                                                   'GCCTCCAGTGAGTCCATACCCCAGGTTTCCTGGTGAGTTCCAAGCCTTTGGAAGCCAGATCTGTG' +
                                                   'ATCCCAAGCCGCCTCCTCCACCGACTGTACTTCATCAACCTTCTCTCATGTTTTTCCGACACTCC' +
                                                   'TGGGTCAGACTCCTNGGTTCATGACTTACTTGCTAGCGTCCCTTCATTTTCCCACAAGTTTGNCC' +
                                                   'CCCACCTCCACTTACCTGTCTGCCCTCAGCTTCTTCCTGGGAGAGAGCCCCCCTTTCACGTAGAC' +
                                                   'ACACCTGGCTGCCTTCTTCACGCCCTGAGGACACTTCTTGGAGATTTCAGAC')
        self.test_target_region.overhang = 250
        self.test_target_region.input_sequence = self.test_input_seq
        self.test_input_seq.strand = "+"
        self.test_target_region.set_primers(min_product_size=300, max_product_size=750, primer_opt_size=20,
                                            primer_min_size=18, primer_max_size=27, primer_opt_tm=60, primer_min_tm=57,
                                            primer_max_tm=63, primer_min_gc=20, primer_max_gc=80)
        # 5 sets of primers should be created
        self.assertEqual(len(self.test_target_region.primers), 5, msg = "Expected 5 primers")

    def test_primer_set_snps(self):
        self.test_target_region.snps = self.test_primer_snp_list
        self.test_primer.set_snps(self.test_target_region)
        forward_snp_ids = []
        for snp in self.test_primer.forward_snps:
            forward_snp_ids.append(snp.snp_id)
        self.assertIn("rsTESTsnp5", forward_snp_ids)
        self.assertIn("rsTESTsnp6", forward_snp_ids)
        self.assertIn("rsTESTsnp7", forward_snp_ids)
        self.assertEqual(len(forward_snp_ids), 3)
        reverse_snp_ids = []
        for snp in self.test_primer.reverse_snps:
            reverse_snp_ids.append(snp.snp_id)
        self.assertIn("rsTESTsnp10", reverse_snp_ids)
        self.assertIn("rsTESTsnp11", reverse_snp_ids)
        self.assertIn("rsTESTsnp12", reverse_snp_ids)
        self.assertEqual(len(reverse_snp_ids), 3)


