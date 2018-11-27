from django.conf.urls import url, include
from django.contrib import admin
from django.views.generic import TemplateView

from autoprimer import views


urlpatterns = [
    url(
        regex = r'^$',
        view = views.home_page,
        name = 'home'),
    url(
        regex = r'^results/$',
        view = TemplateView.as_view(template_name='results.html'),
        name = 'results'),
    url(
        regex = r'^single_target_results/$',
        view = TemplateView.as_view(template_name='single_target_results.html'),
        name = 'single_target_results'),
    url(
        regex = r'^bespoke_target_results/$',
        view = TemplateView.as_view(template_name='bespoke_target_results.html'),
        name = 'bespoke_target_results'),
    url(
        regex = r'^notfound/', 
        view = views.get_gene_info,
        name='notfound.html'),  
    url(
        regex = r'^admin/',
        view = admin.site.urls),
    url(
        regex = r'.csv',
        view = views.download_csv),
    url(
        regex = r'.bed',
        view = views.download_bed),
    url(
        regex = r'^design-by-symbol/',
        view = views.get_parameters_for_symbol,
        name = 'design-by-symbol'),
    url(
        regex = r'^design-by-coordinate/',
        view = views.get_parameters_for_coordinate,
        name = 'design-by-coordinate'),
    url(
        regex = r'^design-for-nipd-bespoke-targets/',
        view = views.get_bespoke_parameters,
        name = 'design-for-nipd-bespoke-targets'),

]