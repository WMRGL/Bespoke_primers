import django_tables2 as tables

class ResultsTable(tables.Table):
    Exon_Target = tables.Column()
    Product_Size = tables.Column()
    Forward_Sequence = tables.Column()
    Forward_Start = tables.Column()
    Forward_Stop = tables.Column()
    Forward_Temp = tables.Column()
    Forward_GC = tables.Column()
    Forward_SNPs = tables.Column()
    Reverse_Sequence = tables.Column()
    Reverse_Start = tables.Column()
    Reverse_Stop = tables.Column()
    Reverse_Temp = tables.Column()
    Reverse_GC = tables.Column()
    Reverse_SNPs = tables.Column()

    class Meta:
        attrs = {'class': 'table table-bordered'}
        orderable = False

    def make_table(input_seq):
        data = []
        for target in input_seq.target_regions:
            for primer in target.primers:
                forward_primer_snps = []
                reverse_primer_snps = []
                for snp in primer.forward_snps:
                    forward_primer_snps.append(snp.snp_id + ' (' + str(round(snp.av_het,4)) + ')\n')
                f_snps = ''.join(forward_primer_snps)
                for snp in primer.reverse_snps:
                    reverse_primer_snps.append(snp.snp_id + ' (' + str(round(snp.av_het,4)) + ')\n')
                r_snps = ''.join(reverse_primer_snps)
                data.append({'Exon_Target': target.target_id,
                             'Product_Size': primer.product_size,
                             'Forward_Sequence': primer.forward_seq,
                             'Forward_Start': primer.forward_genomic_coords[0],
                             'Forward_Stop': primer.forward_genomic_coords[1],
                             'Forward_Temp': round(primer.forward_tm, 2),
                             'Forward_GC': round(primer.forward_gc, 2),
                             'Forward_SNPs': f_snps,
                             'Reverse_Sequence': primer.reverse_seq,
                             'Reverse_Start': primer.reverse_genomic_coords[0],
                             'Reverse_Stop': primer.reverse_genomic_coords[1],
                             'Reverse_Temp': round(primer.reverse_tm, 2),
                             'Reverse_GC': round(primer.reverse_gc, 2),
                             'Reverse_SNPs': r_snps,
                             })
        table = ResultsTable(data)
        return table

class SingleTargetResultsTable(tables.Table):
    Product_Size = tables.Column()
    Forward_Sequence = tables.Column()
    Forward_Start = tables.Column()
    Forward_Stop = tables.Column()
    Forward_Temp = tables.Column()
    Forward_GC = tables.Column()
    Forward_SNPs = tables.Column()
    Reverse_Sequence = tables.Column()
    Reverse_Start = tables.Column()
    Reverse_Stop = tables.Column()
    Reverse_Temp = tables.Column()
    Reverse_GC = tables.Column()
    Reverse_SNPs = tables.Column()

    class Meta:
        attrs = {'class': 'table table-bordered'}
        orderable = False

    def make_table(input_seq):
        data = []
        for target in input_seq.target_regions:
            for primer in target.primers:
                forward_primer_snps = []
                reverse_primer_snps = []
                for snp in primer.forward_snps:
                    forward_primer_snps.append(snp.snp_id + ' (' + str(round(snp.av_het,4)) + ')\n')
                f_snps = ''.join(forward_primer_snps)
                for snp in primer.reverse_snps:
                    reverse_primer_snps.append(snp.snp_id + ' (' + str(round(snp.av_het,4)) + ')\n')
                r_snps = ''.join(reverse_primer_snps)
                data.append({'Product_Size': primer.product_size,
                             'Forward_Sequence': primer.forward_seq,
                             'Forward_Start': primer.forward_genomic_coords[0],
                             'Forward_Stop': primer.forward_genomic_coords[1],
                             'Forward_Temp': round(primer.forward_tm, 2),
                             'Forward_GC': round(primer.forward_gc, 2),
                             'Forward_SNPs': f_snps,
                             'Reverse_Sequence': primer.reverse_seq,
                             'Reverse_Start': primer.reverse_genomic_coords[0],
                             'Reverse_Stop': primer.reverse_genomic_coords[1],
                             'Reverse_Temp': round(primer.reverse_tm, 2),
                             'Reverse_GC': round(primer.reverse_gc, 2),
                             'Reverse_SNPs': r_snps,
                             })
        table = SingleTargetResultsTable(data)
        return table


class BespokeTargetResultsTable(tables.Table):
    Product_Size = tables.Column()
    Forward_Sequence = tables.Column()
    Forward_Start = tables.Column()
    Forward_Stop = tables.Column()
    Forward_Temp = tables.Column()
    Forward_GC = tables.Column()
    Forward_SNPs = tables.Column()
    Reverse_Sequence = tables.Column()
    Reverse_Start = tables.Column()
    Reverse_Stop = tables.Column()
    Reverse_Temp = tables.Column()
    Reverse_GC = tables.Column()
    Reverse_SNPs = tables.Column()

    class Meta:
        attrs = {'class': 'table table-bordered'}
        orderable = False

    def make_table(input_seq):
        data = []
        for target in input_seq.target_regions:
            for primer in target.primers:
                forward_primer_snps = []
                reverse_primer_snps = []
                for snp in primer.forward_snps:
                    forward_primer_snps.append(snp.snp_id + ' (' + str(round(snp.av_het,4)) + ')\n')
                f_snps = ''.join(forward_primer_snps)
                for snp in primer.reverse_snps:
                    reverse_primer_snps.append(snp.snp_id + ' (' + str(round(snp.av_het,4)) + ')\n')
                r_snps = ''.join(reverse_primer_snps)
                data.append({'Product_Size': primer.product_size,
                             'Forward_Sequence': primer.forward_seq,
                             'Forward_Start': primer.forward_genomic_coords[0],
                             'Forward_Stop': primer.forward_genomic_coords[1],
                             'Forward_Temp': round(primer.forward_tm, 2),
                             'Forward_GC': round(primer.forward_gc, 2),
                             'Forward_SNPs': f_snps,
                             'Reverse_Sequence': primer.reverse_seq,
                             'Reverse_Start': primer.reverse_genomic_coords[0],
                             'Reverse_Stop': primer.reverse_genomic_coords[1],
                             'Reverse_Temp': round(primer.reverse_tm, 2),
                             'Reverse_GC': round(primer.reverse_gc, 2),
                             'Reverse_SNPs': r_snps,
                             })
        table = BespokeTargetResultsTable(data)
        return table