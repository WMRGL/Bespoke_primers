from django import forms

GENOME_BUILD_CHOICES=[('GRCh37', 'GRCh37'),('GRCh38', 'GRCh38')]
CHROMOSOME_CHOICES=[
    ('chr1', 'chr1'),
    ('chr2', 'chr2'),
    ('chr3', 'chr3'),
    ('chr4', 'chr4'),
    ('chr5', 'chr5'),
    ('chr6', 'chr6'),
    ('chr7', 'chr7'),
    ('chr8', 'chr8'),
    ('chr9', 'chr9'),
    ('chr10', 'chr10'),
    ('chr11', 'chr11'),
    ('chr12', 'chr12'),
    ('chr13', 'chr13'),
    ('chr14', 'chr14'),
    ('chr15', 'chr15'),
    ('chr16', 'chr16'),
    ('chr17', 'chr17'),
    ('chr18', 'chr18'),
    ('chr19', 'chr19'),
    ('chr20', 'chr20'),
    ('chr21', 'chr21'),
    ('chr22', 'chr22'),
    ('chrX', 'chrX'),
    ('chrY', 'chrY'),
]

class SymbolDesignForm(forms.Form):
    genome_build = forms.CharField(label="Reference Genome",
                                     widget=forms.RadioSelect(choices=GENOME_BUILD_CHOICES))
    gene_symbol = forms.CharField(label="Gene symbol")
    max_snp_avhet = forms.FloatField(label="Maximum SNP avHet", initial=0.02, max_value=0.5, min_value=0)
    min_product_size = forms.IntegerField(label='Minimum product size', initial=300)
    max_product_size = forms.IntegerField(label='Maximum product size', initial=750, min_value=300)
    primer_opt_size = forms.IntegerField(label="Optimum primer size", initial=20)
    primer_min_size = forms.IntegerField(label="Minimum primer size", initial=18)
    primer_max_size = forms.IntegerField(label="Maximum primer size", initial=27)
    primer_opt_tm = forms.FloatField(label="Optimum annealing temp", initial=60)
    primer_min_tm = forms.FloatField(label="Minimum annealing temp", initial=57)
    primer_max_tm = forms.FloatField(label="Maximum annealing temp", initial=63)
    primer_min_gc = forms.FloatField(label="Minimum GC percent", initial=20, max_value=100, min_value=0)
    primer_max_gc = forms.FloatField(label="Maximum GC percent", initial=80, max_value=100, min_value=0)

    def __init__(self, *args, **kwargs):
        super(SymbolDesignForm, self).__init__(*args, **kwargs)
        for field_name, field in self.fields.items():
            field.widget.attrs['class'] = 'form-horizontal'
            field.widget.attrs['size'] = 1

    def clean(self):
        cleaned_data = super(SymbolDesignForm, self).clean()
        max_product_size = cleaned_data.get('max_product_size')
        min_product_size = cleaned_data.get('min_product_size')
        primer_opt_size = cleaned_data.get('primer_opt_size')
        primer_min_size = cleaned_data.get('primer_min_size')
        primer_max_size = cleaned_data.get('primer_max_size')
        primer_opt_tm = cleaned_data.get('primer_opt_tm')
        primer_min_tm = cleaned_data.get('primer_min_tm')
        primer_max_tm = cleaned_data.get('primer_max_tm')
        primer_min_gc = cleaned_data.get('primer_min_gc')
        primer_max_gc = cleaned_data.get('primer_max_gc')
        if max_product_size < min_product_size:
            self.add_error('max_product_size',
                           "Maximum product size must be greater than or equal to minimum product size")
            self.add_error('min_product_size',
                           "Minimum product size must be less than or equal to maximum product size")
        if primer_max_size < primer_min_size:
            self.add_error('primer_max_size',
                           "Maximum primer size must be greater than or equal to minimum primer size")
            self.add_error('primer_min_size',
                           "Minimum primer size must be less than or equal to maximum primer size")
        if primer_min_size > primer_opt_size or primer_max_size < primer_opt_size:
            self.add_error('primer_opt_size',
                           ("Optimum primer size must be greater than or equal to minimum primer size and less "
                            "than or equal to maximum primer size"))
        if min_product_size < primer_max_size:
            self.add_error('min_product_size', 'Minimum product size must be greater than maximum primer size')
        if primer_max_tm < primer_min_tm:
            self.add_error('primer_max_tm',
                           ("Maximum annealing temperature must be greater than or equal to minimum annealing "
                            "temperature"))
            self.add_error('primer_min_tm',
                           "Minimum annealing temperature must be less than or equal to maximum annealing temperature")
        if primer_min_tm > primer_opt_tm or primer_max_tm < primer_opt_tm:
            self.add_error('primer_opt_tm',
                           ("Optimum annealing temperature must be greater than or equal to minimum annealing "
                            "temperature and less than or equal to maximum annealing temperature"))
        if primer_max_gc < primer_min_gc:
            self.add_error('primer_max_gc',
                           "Maximum GC percent must be greater than or equal to minimum GC percent")
            self.add_error('primer_min_gc',
                           "Minimum GC percent must be less than or equal to maximum GC percent")

class SingleTargetDesignForm(forms.Form):
    genome_build = forms.CharField(label="Reference Genome",
                                     widget=forms.RadioSelect(choices=GENOME_BUILD_CHOICES))
    chromosome = forms.ChoiceField(label="Chromosome", choices=CHROMOSOME_CHOICES)
    coordinate = forms.IntegerField(label='Coordinate')
    max_snp_avhet = forms.FloatField(label="Maximum SNP avHet", initial=0.02, max_value=0.5, min_value=0)
    min_product_size = forms.IntegerField(label='Minimum product size', initial=300)
    max_product_size = forms.IntegerField(label='Maximum product size', initial=750, min_value=300)
    primer_opt_size = forms.IntegerField(label="Optimum primer size", initial=20)
    primer_min_size = forms.IntegerField(label="Minimum primer size", initial=18)
    primer_max_size = forms.IntegerField(label="Maximum primer size", initial=27)
    primer_opt_tm = forms.FloatField(label="Optimum annealing temp", initial=60)
    primer_min_tm = forms.FloatField(label="Minimum annealing temp", initial=57)
    primer_max_tm = forms.FloatField(label="Maximum annealing temp", initial=63)
    primer_min_gc = forms.FloatField(label="Minimum GC percent", initial=20, max_value=100, min_value=0)
    primer_max_gc = forms.FloatField(label="Maximum GC percent", initial=80, max_value=100, min_value=0)

    def __init__(self, *args, **kwargs):
        super(SingleTargetDesignForm, self).__init__(*args, **kwargs)
        for field_name, field in self.fields.items():
            field.widget.attrs['class'] = 'form-horizontal'
            field.widget.attrs['size'] = 1

    def clean(self):
        cleaned_data = super(SingleTargetDesignForm, self).clean()
        coordinate = cleaned_data.get('coordinate')
        max_product_size = cleaned_data.get('max_product_size')
        min_product_size = cleaned_data.get('min_product_size')
        primer_opt_size = cleaned_data.get('primer_opt_size')
        primer_min_size = cleaned_data.get('primer_min_size')
        primer_max_size = cleaned_data.get('primer_max_size')
        primer_opt_tm = cleaned_data.get('primer_opt_tm')
        primer_min_tm = cleaned_data.get('primer_min_tm')
        primer_max_tm = cleaned_data.get('primer_max_tm')
        primer_min_gc = cleaned_data.get('primer_min_gc')
        primer_max_gc = cleaned_data.get('primer_max_gc')
        if coordinate < 0:
            self.add_error('coordinate',
                           "Coordinate must be a postive integer")
        if max_product_size < min_product_size:
            self.add_error('max_product_size',
                           "Maximum product size must be greater than or equal to minimum product size")
            self.add_error('min_product_size',
                           "Minimum product size must be less than or equal to maximum product size")
        if primer_max_size < primer_min_size:
            self.add_error('primer_max_size',
                           "Maximum primer size must be greater than or equal to minimum primer size")
            self.add_error('primer_min_size',
                           "Minimum primer size must be less than or equal to maximum primer size")
        if primer_min_size > primer_opt_size or primer_max_size < primer_opt_size:
            self.add_error('primer_opt_size',
                           ("Optimum primer size must be greater than or equal to minimum primer size and less "
                            "than or equal to maximum primer size"))
        if min_product_size < primer_max_size:
            self.add_error('min_product_size', 'Minimum product size must be greater than maximum primer size')
        if primer_max_tm < primer_min_tm:
            self.add_error('primer_max_tm',
                           ("Maximum annealing temperature must be greater than or equal to minimum annealing "
                            "temperature"))
            self.add_error('primer_min_tm',
                           "Minimum annealing temperature must be less than or equal to maximum annealing temperature")
        if primer_min_tm > primer_opt_tm or primer_max_tm < primer_opt_tm:
            self.add_error('primer_opt_tm',
                           ("Optimum annealing temperature must be greater than or equal to minimum annealing "
                            "temperature and less than or equal to maximum annealing temperature"))
        if primer_max_gc < primer_min_gc:
            self.add_error('primer_max_gc',
                           "Maximum GC percent must be greater than or equal to minimum GC percent")
            self.add_error('primer_min_gc',
                           "Minimum GC percent must be less than or equal to maximum GC percent")


class BespokeTargetDesignForm(forms.Form):
    genome_build = forms.CharField(label="Reference Genome",
                                     widget=forms.RadioSelect(choices=GENOME_BUILD_CHOICES))
    chromosome = forms.ChoiceField(label="Chromosome", choices=CHROMOSOME_CHOICES)
    coordinate = forms.IntegerField(label='Coordinate')
    # add an end coordinate
    max_snp_avhet = forms.FloatField(label="Maximum SNP avHet", initial=0.02, max_value=0.5, min_value=0)
    min_product_size = forms.IntegerField(label='Minimum product size', initial=50)
    max_product_size = forms.IntegerField(label='Maximum product size', initial=90, min_value=50)
    primer_opt_size = forms.IntegerField(label="Optimum primer size", initial=22)
    primer_min_size = forms.IntegerField(label="Minimum primer size", initial=18)
    primer_max_size = forms.IntegerField(label="Maximum primer size", initial=30)
    primer_opt_tm = forms.FloatField(label="Optimum annealing temp", initial=62)
    primer_min_tm = forms.FloatField(label="Minimum annealing temp", initial=60)
    primer_max_tm = forms.FloatField(label="Maximum annealing temp", initial=65)
    primer_min_gc = forms.FloatField(label="Minimum GC percent", initial=20, max_value=100, min_value=0)
    primer_max_gc = forms.FloatField(label="Maximum GC percent", initial=80, max_value=100, min_value=0)

    def __init__(self, *args, **kwargs):
        super(BespokeTargetDesignForm, self).__init__(*args, **kwargs)
        for field_name, field in self.fields.items():
            field.widget.attrs['class'] = 'form-horizontal'
            field.widget.attrs['size'] = 1

    def clean(self):
        cleaned_data = super(BespokeTargetDesignForm, self).clean()
        coordinate = cleaned_data.get('coordinate')
        max_product_size = cleaned_data.get('max_product_size')
        min_product_size = cleaned_data.get('min_product_size')
        primer_opt_size = cleaned_data.get('primer_opt_size')
        primer_min_size = cleaned_data.get('primer_min_size')
        primer_max_size = cleaned_data.get('primer_max_size')
        primer_opt_tm = cleaned_data.get('primer_opt_tm')
        primer_min_tm = cleaned_data.get('primer_min_tm')
        primer_max_tm = cleaned_data.get('primer_max_tm')
        primer_min_gc = cleaned_data.get('primer_min_gc')
        primer_max_gc = cleaned_data.get('primer_max_gc')
        if coordinate < 0:
            self.add_error('coordinate',
                           "Coordinate must be a postive integer")
        if max_product_size < min_product_size:
            self.add_error('max_product_size',
                           "Maximum product size must be greater than or equal to minimum product size")
            self.add_error('min_product_size',
                           "Minimum product size must be less than or equal to maximum product size")
        if primer_max_size < primer_min_size:
            self.add_error('primer_max_size',
                           "Maximum primer size must be greater than or equal to minimum primer size")
            self.add_error('primer_min_size',
                           "Minimum primer size must be less than or equal to maximum primer size")
        if primer_min_size > primer_opt_size or primer_max_size < primer_opt_size:
            self.add_error('primer_opt_size',
                           ("Optimum primer size must be greater than or equal to minimum primer size and less "
                            "than or equal to maximum primer size"))
        if min_product_size < primer_max_size:
            self.add_error('min_product_size', 'Minimum product size must be greater than maximum primer size')
        if primer_max_tm < primer_min_tm:
            self.add_error('primer_max_tm',
                           ("Maximum annealing temperature must be greater than or equal to minimum annealing "
                            "temperature"))
            self.add_error('primer_min_tm',
                           "Minimum annealing temperature must be less than or equal to maximum annealing temperature")
        if primer_min_tm > primer_opt_tm or primer_max_tm < primer_opt_tm:
            self.add_error('primer_opt_tm',
                           ("Optimum annealing temperature must be greater than or equal to minimum annealing "
                            "temperature and less than or equal to maximum annealing temperature"))
        if primer_max_gc < primer_min_gc:
            self.add_error('primer_max_gc',
                           "Maximum GC percent must be greater than or equal to minimum GC percent")
            self.add_error('primer_min_gc',
                           "Minimum GC percent must be less than or equal to maximum GC percent")

