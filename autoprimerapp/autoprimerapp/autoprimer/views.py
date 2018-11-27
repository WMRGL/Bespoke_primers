import csv
import os
from django.shortcuts import render, redirect, render_to_response
from django.template import RequestContext

from django.http import HttpResponseRedirect, HttpResponse

from .forms import SymbolDesignForm, SingleTargetDesignForm, BespokeTargetDesignForm
from .tables import ResultsTable, SingleTargetResultsTable, BespokeTargetResultsTable

from autoprimer.autoprimer import design_by_symbol
from autoprimer.singletarget import design_by_coordinate
from autoprimer.bespoketarget import bespoke_design

def home_page(request):
    return render(request, 'index.html')

def get_parameters_for_symbol(request):
    if request.method == 'POST':
        form = SymbolDesignForm(request.POST)
        if form.is_valid():
            genome_build = form.cleaned_data.get("genome_build")
            gene_symbol = form.cleaned_data.get("gene_symbol").upper()
            gene_info = get_gene_info(gene_symbol, genome_build)
            if gene_info == None:
                return render(request, 'notfound.html', {'gene_symbol': gene_symbol})
            ng_number = gene_info[0]
            chromosome = gene_info[1]
            start = int(gene_info[2])
            strand = gene_info[4]
            print(genome_build, start)
            max_avhet = form.cleaned_data.get("max_snp_avhet")
            min_product_size = form.cleaned_data.get("min_product_size")
            max_product_size = form.cleaned_data.get("max_product_size")
            primer_opt_size = form.cleaned_data.get("primer_opt_size")
            primer_min_size = form.cleaned_data.get("primer_min_size")
            primer_max_size = form.cleaned_data.get("primer_max_size")
            primer_opt_tm = form.cleaned_data.get("primer_opt_tm")
            primer_min_tm = form.cleaned_data.get("primer_min_tm")
            primer_max_tm = form.cleaned_data.get("primer_max_tm")
            primer_min_gc = form.cleaned_data.get("primer_min_gc")
            primer_max_gc = form.cleaned_data.get("primer_max_gc")
            input_seq = design_by_symbol(genome_build, ng_number, chromosome, start, strand, max_avhet,
                                         min_product_size, max_product_size, primer_opt_size, primer_min_size,
                                         primer_max_size, primer_opt_tm, primer_min_tm, primer_max_tm, primer_min_gc,
                                         primer_max_gc)
            refseq_id = input_seq.refseq_id
            gene = input_seq.gene_name
            chromosome = input_seq.chrom_number
            refseq_coord_start = input_seq.genomic_coords[0]
            refseq_coord_stop = input_seq.genomic_coords[1]
            results = ResultsTable.make_table(input_seq)
            return render(request, 'results.html', {
                'genome_build': genome_build,
                'refseq_id': refseq_id,
                'gene': gene,
                'chromosome': chromosome,
                'refseq_coord_start': refseq_coord_start,
                'refseq_coord_stop': refseq_coord_stop,
                'results': results,
            })
    else:
        form = SymbolDesignForm()
    return render(request, 'design-by-symbol.html', {'form': form})

def get_parameters_for_coordinate(request):
    if request.method == 'POST':
        form = SingleTargetDesignForm(request.POST)
        if form.is_valid():
            genome_build = form.cleaned_data.get("genome_build")
            chromosome = form.cleaned_data.get("chromosome")
            coordinate = form.cleaned_data.get("coordinate")
            max_avhet = form.cleaned_data.get("max_snp_avhet")
            min_product_size = form.cleaned_data.get("min_product_size")
            max_product_size = form.cleaned_data.get("max_product_size")
            primer_opt_size = form.cleaned_data.get("primer_opt_size")
            primer_min_size = form.cleaned_data.get("primer_min_size")
            primer_max_size = form.cleaned_data.get("primer_max_size")
            primer_opt_tm = form.cleaned_data.get("primer_opt_tm")
            primer_min_tm = form.cleaned_data.get("primer_min_tm")
            primer_max_tm = form.cleaned_data.get("primer_max_tm")
            primer_min_gc = form.cleaned_data.get("primer_min_gc")
            primer_max_gc = form.cleaned_data.get("primer_max_gc")
            input_seq = design_by_coordinate(genome_build, chromosome, coordinate, max_avhet, min_product_size,
                                             max_product_size, primer_opt_size, primer_min_size, primer_max_size,
                                             primer_opt_tm, primer_min_tm, primer_max_tm, primer_min_gc,
                                             primer_max_gc)
            chromosome = input_seq.chrom_number
            results = SingleTargetResultsTable.make_table(input_seq)
            filename = chromosome + '-' + str(coordinate)
            return render(request, 'single_target_results.html', {
                'filename': filename,
                'genome_build': genome_build,
                'chromosome': chromosome,
                'coordinate': coordinate,
                'results': results,
            })
    else:
        form = SingleTargetDesignForm()
    return render(request, 'design-by-coordinate.html', {'form': form})

def get_bespoke_parameters(request):
    if request.method == 'POST':
        form = BespokeTargetDesignForm(request.POST)
        if form.is_valid():
            genome_build = form.cleaned_data.get("genome_build")
            chromosome = form.cleaned_data.get("chromosome")
            coordinate = form.cleaned_data.get("coordinate")
            max_avhet = form.cleaned_data.get("max_snp_avhet")
            min_product_size = form.cleaned_data.get("min_product_size")
            max_product_size = form.cleaned_data.get("max_product_size")
            primer_opt_size = form.cleaned_data.get("primer_opt_size")
            primer_min_size = form.cleaned_data.get("primer_min_size")
            primer_max_size = form.cleaned_data.get("primer_max_size")
            primer_opt_tm = form.cleaned_data.get("primer_opt_tm")
            primer_min_tm = form.cleaned_data.get("primer_min_tm")
            primer_max_tm = form.cleaned_data.get("primer_max_tm")
            primer_min_gc = form.cleaned_data.get("primer_min_gc")
            primer_max_gc = form.cleaned_data.get("primer_max_gc")
            input_seq = bespoke_design(genome_build, chromosome, coordinate, max_avhet, min_product_size,
                                             max_product_size, primer_opt_size, primer_min_size, primer_max_size,
                                             primer_opt_tm, primer_min_tm, primer_max_tm, primer_min_gc,
                                             primer_max_gc)
            chromosome = input_seq.chrom_number
            results = BespokeTargetResultsTable.make_table(input_seq)
            filename = chromosome + '-' + str(coordinate)
            return render(request, 'bespoke_target_results.html', {
                'filename': filename,
                'genome_build': genome_build,
                'chromosome': chromosome,
                'coordinate': coordinate,
                'results': results,
            })
    else:
        form = BespokeTargetDesignForm()
    return render(request, 'design-for-nipd-bespoke-targets.html', {'form': form})

def get_gene_info(gene_symbol, genome_build):
    ref_seq_dict = {}

    if genome_build == "GRCh37":
        with open("input/RefSeqGenes37.txt") as tsv:
            for line in csv.reader(tsv, delimiter=" "):
                ref_seq_dict[line[0]] = [line[1],line[2],line[3],line[4],line[5]]
    elif genome_build == "GRCh38":
        with open("input/RefSeqGenes38.txt") as tsv:
            for line in csv.reader(tsv, delimiter=" "):
                ref_seq_dict[line[0]] = [line[1],line[2],line[3],line[4],line[5]]

    result = None
    if gene_symbol in ref_seq_dict:
        result = ref_seq_dict[gene_symbol]
    return result


def download_csv(request):
    filename = request.path.split('/')[-1]
    # with open("/srv/primer_design/s_drive/designs/" + filename, 'rb') as fh:
    with open("/media/sf_S_DRIVE/genomic_resources/primer_design/designs/" + filename, 'rb') as fh:
        response = HttpResponse(fh.read(), content_type="application/vnd.ms-excel")
        response['Content-Disposition'] = 'inline; filename=' + os.path.basename('output/' + filename)
        return response

def download_bed(request):
    filename = request.path.split('/')[-1]
    # with open("/srv/primer_design/s_drive/designs/" + filename, 'rb') as fh:
    with open("/media/sf_S_DRIVE/genomic_resources/primer_design/designs/" + filename, 'rb') as fh:
        response = HttpResponse(fh.read(), content_type="application/vnd.ms-excel")
        response['Content-Disposition'] = 'inline; filename=' + os.path.basename('output/' + filename)
        return response
