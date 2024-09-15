from django import forms
from django.shortcuts import render
from base import forms
from base.curves import *
from sympy import nextprime
import math
from gmpy2 import mpz
# Create your views here.
# a = 0
# d = 0
# p = 0
# new_p = 0
# set = False

def home(request):
    # global a,d,p,new_p,set
    new_p = 0
    prime = 0
    
    adp_form = forms.adp_form()

    if request.method == "POST":

        adp_form = forms.adp_form(request.POST)

        if adp_form.is_valid():
            
            opt = request.session['opt'] = adp_form.cleaned_data['opt']
            a = request.session['a'] = adp_form.cleaned_data['a']
            d = request.session['d'] = adp_form.cleaned_data['d']
            p = request.session['p'] = adp_form.cleaned_data['p']
            request.session['set'] = True
            new_p = request.session['new_p'] = nextprime(p-1)
            lo = request.session['lo'] = int(new_p + 1 - 2*(new_p**0.5))
            hi = request.session['hi'] = int(new_p + 1 + 2*(new_p**0.5))
            # lo = request.session['lo'] = mpz(new_p + 1 - 2*(new_p**0.5))
            # hi = request.session['hi'] = mpz(new_p + 1 + 2*(new_p**0.5))
            prime = (new_p == p)

            #deciding on labels        
            a_label = 'a'
            d_label = 'd'
            p_label = 'p'

            if opt == '3':              
              a_label = 'A'
              d_label = 'B'
                        
            return render(request, 'base/home.html', {'adp_form': adp_form, 'stage': 2, 'a': a, 'd': d, 'p': p, 'new_p': new_p, 'lo': lo, 'hi': hi, 'prime': prime, 'a_label': a_label, 'd_label': d_label, 'p_label': p_label})
    return render(request, 'base/home.html', {'adp_form': adp_form, 'stage': 1})

def calc(request, start=0):
    # handle start value
    start = int(start)

    # global a,d,p,new_p,set
    if not request.session['set']:
        return render(request,'base/notset.html')
    else:
        curve = t_edwards
        opt1 = request.session['opt']
        # print("views -> opt1", opt1)

        a = request.session['a']
        d = request.session['d']
        new_p = request.session['new_p']
        curve = t_edwards
        if opt1 == '1':
            curve = t_edwards
        elif opt1 == '2':
            curve = s_weirstrass_curve
        elif opt1 == '3':
            curve = montgomery_curve
        points = curve.generatePoints(a, d, new_p, start)
        opt_form = forms.opt_form()

        a_label = 'a'
        d_label = 'd'
        p_label = 'p'
        
        if opt1 == '3':              
          a_label = 'A'
          d_label = 'B'
        
        # Operations +,-,* trigger POST
        if request.method == "POST":

            opt_form = forms.opt_form(request.POST)
            
            if opt_form.is_valid():

                opt = opt_form.cleaned_data['opt']
                x1 = opt_form.cleaned_data['x1']
                y1 = opt_form.cleaned_data['y1']
                x2 = opt_form.cleaned_data['x2']
                y2 = opt_form.cleaned_data['y2']

                # print(x1,y1,x2,y2)

                x_res = 0
                y_res = 0
                k = 0

                if(opt == '2'):
                    (x_res,y_res) = curve.addpoints(a,d,new_p,(x1,y1), (x2,y2))
                elif(opt == '3'):
                    (x_res,y_res) = curve.substractpoints(a,d,new_p,(x1,y1), (x2,y2))
                elif(opt == '4'):
                    (x_res,y_res) = curve.doublepoint(a,d,new_p,(x1,y1))
                elif(opt == '5'):
                    (x_res,y_res) = curve.multiplypoint(a,d,new_p,(x1,y1), x2)
                elif(opt == '6'):
                    # (x_res, y_res) = curve.bsgs(a,d,new_p,(x1,y1),(x2,y2))
                    k = curve.bsgs(a,d,new_p,(x1,y1),(x2,y2))

                return render(request,'base/calculate.html',{'opt_form': opt_form, 'a': a, 'd': d, 'p': new_p, 'xarray': points[0], 'yarray': points[1], 'Array': zip(points[0], points[1]), 'point_count': len(points[0]), 'x_res': x_res, 'y_res': y_res, 'k':k, 'result': True, 'start': start, 'end': min(new_p-1, start+999), 'prev': max(0, start-1000), 'next': min(new_p-1, start+1000), 'p_minus_1': new_p-1,'curve': opt1, 'a_label': a_label, 'd_label': d_label, 'p_label': p_label})

        # GET
        return render(request,'base/calculate.html',{'opt_form': opt_form, 'a': a, 'd': d, 'p': new_p, 'xarray': points[0], 'yarray': points[1], 'Array': zip(points[0], points[1]), 'point_count': len(points[0]), 'start': start, 'end': min(new_p-1, start+999), 'prev': max(0, start-1000), 'next': min(new_p-1, start+1000), 'p_minus_1': new_p-1,'curve': opt1, 'a_label': a_label, 'd_label': d_label, 'p_label': p_label})
    
def credits(request):
    return render(request, 'base/credits.html')