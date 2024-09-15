from django.urls import path
from . import views

urlpatterns = [
    path('',views.home,name="home"),
    path('calculate/<str:start>/',views.calc,name="calculate"),
    path('credits/', views.credits, name="credits")
]