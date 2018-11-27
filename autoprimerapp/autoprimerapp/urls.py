
from django.conf.urls import url, include
from django.contrib import admin
from django.views.generic import TemplateView

from autoprimer import views


urlpatterns = [
    url(r'^$',
        TemplateView.as_view(template_name='index.html'),
        name='home'),
    url(r'^results/$',
        TemplateView.as_view(template_name='results.html'),
        name='results'),
    url(r'^single_target_results/$',
        TemplateView.as_view(template_name='single_target_results.html'),
        name='single_target_results'),
    url(r'^bespoke_target_results/$',
        TemplateView.as_view(template_name='bespoke_target_results.html'),
        name='bespoke_target_results'),
    url(r'^notfound/', views.get_gene_info, name='notfound.html'),
    url(r'^admin/', admin.site.urls),
    url(r'.csv', views.download_csv),
    url(r'.bed', views.download_bed),
    url(r'^design-by-symbol/', views.get_parameters_for_symbol, name='design-by-symbol'),
    url(r'^design-by-coordinate/', views.get_parameters_for_coordinate, name='design-by-coordinate'),

]

