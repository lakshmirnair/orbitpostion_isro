from __future__ import print_function, division
from PyAstronomy import pyasl
from werkzeug.utils import secure_filename

from win10toast import ToastNotifier

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from statsmodels.tsa.arima_model import ARIMA
import math
import os
import flask
from flask import Flask, render_template, request

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = '/static/plot'


def add(file):
    data = []
    f = open(file, 'r')
    lines = f.readlines()
    k = 7
    n = 7
    i = -1
    row = []
    sss = file[:3]
    date = file[3:5] + '/' + file[5:7] + '/20' + file[7:9]
    row.append(sss)
    row.append(date)
    for line in lines[7:]:
        temp = line.split(' ')
        temp = list(filter(None, temp))
        temp = [w.replace('D', 'e') for w in temp]
        if k == n:
            n = n + 8
            i = i + 1
            trow = row.copy()
            t = temp.copy()
            t[-1] = t[-1].rstrip()
            trow.extend(t)
            data.append(trow)
        elif k > 7:
            t = temp.copy()
            t[-1] = t[-1].rstrip()
            data[i].extend(t)
        k = k + 1
    return data

@app.route('/')
def hello():
    for col in df.columns:
        filename = str(col)+'.png'
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        # full_filename="static/plot/pic.png"
        html='<img src="hello.png" />'
        return render_template("root.html", image=full_filename)


@app.route('/front', methods=['GET', 'POST'])
def front():
    if request.method == 'GET':
        return render_template("front.html")
    else:
        ks = pyasl.MarkleyKESolver()

        sat_no = request.form['sat_no']
        df = pd.read_csv('sornew.csv')
        df = df.loc[df['PRN'] == int(sat_no)]
        df = df[['Del_n', 'e', 'sqrt(A)', 'i0', 'omega', 'M0']]
        M0 = df['M0'].values[0]
        e = df['e'].values[0]
        print("Eccentric anomaly: ", ks.getE(M0, e))



        for col in df.columns:
            df1 = df[[col]]
            # df1 = df1.loc[:,len(df1)]
            plt.xlabel(col)
            plt.plot(df1)
            plt.savefig('static/plot/' + col + '.png', format='png')
            plt.cla()

        Image_folder = os.path.join('static', 'plot')

        return render_template("root.html")


toaster = ToastNotifier()


@app.route('/uploader', methods=['POST'])
def uploader():
    mdata = []
    k=0
    #f = request.files['file']
    uploaded_files = flask.request.files.getlist("file")
    for j in uploaded_files:
        file = j.filename
        j.save(secure_filename(j.filename))
        data = add(j.filename)
        columns = ['SSS', 'Date',
                   'PRN', 'Epoch-year', 'Epoch-Month', 'Epoch-Day', 'Epoch-Hour', 'Epoch-Minute', 'Epoch-Sec',
                   'SV_Clock_Bias', 'SV_Clock_Drift', 'SV_Clock_Drift_Rate',
               'IODE', 'Crs', 'Del_n', 'M0', 'Cuc', 'e', 'Cus', 'sqrt(A)', 'Toe', 'Cic', 'OMEGA', 'Cis',
               'i0', 'Crc', 'omega', 'OMEGA_dot', 'I_dot', 'Codes', 'GPS_week', 'L2_P_Data_flag',
               'SV_accuracy', 'SV_health', 'Tgd', 'IODC', 'T_Tx', 'Fit_Interval', '']
        df = pd.DataFrame(data, columns=columns)
        df.sort_values(['PRN', 'Epoch-year', 'Epoch-Month', 'Epoch-Day', 'Epoch-Hour', 'Epoch-Minute', 'Epoch-Sec'],
                 ascending=[True, True, True, True, True, True, True], inplace=True)
        df.drop(df.columns[len(df.columns) - 1], axis=1, inplace=True)
    #vdf = pd.read_csv('sornew.csv')
    #vdf = vdf.append(df, ignore_index=True)
    #vdf.sort_values(['PRN', 'Epoch-year', 'Epoch-Month', 'Epoch-Day', 'Epoch-Hour', 'Epoch-Minute', 'Epoch-Sec'],
         #           ascending=[True, True, True, True, True, True, True], inplace=True)
        df.to_csv('sornew1.csv', index=False,mode='a')
        # vdf.to_csv('sornew1.csv', mode='a', header=False)
    #vdf.to_csv('sornew1.csv', mode='w')

    toaster.show_toast("File uploaded successfully", "PLEASE SELECT SATELLITE NUMBER")
    return render_template('front.html')


@app.route('/mod', methods=['GET', 'POST'])
def mod():
    if request.method == 'GET':
        return render_template("mod.html")
    else:
        sat_no = request.form['sat_no']
        fromhistory = request.form['fromfiles_date']
        tohistory = request.form['tofiles_date']
        toreq = request.form['toreq_date']
        #timestamp = int(request.form['timestamp'])
        fromhistory = int(fromhistory[-2:])
        tohistory = int(tohistory[-2:])
        toreq=int(toreq[-2:])
        fromreq= tohistory+1
        t = toreq-fromreq
        timestamp = (t+1) * 20
        print(timestamp)
        m=[]
        order = [[0,0,1],[5,2,0],[0,0,1],[5,1,1],[5,2,1],[5,1,0],[5, 1, 0], [6, 1, 0], [4, 1, 0], [5, 1, 0],[5,0,0],[0,1,1]]
        cols = ['Cus','Crc','Cis','Cic','Cuc','Crs','M0', 'OMEGA', 'e', 'Del_n','Toe','omega']
        head=['SSN','PRN','Year','Month','Day','Hour','Minutes','Seconds','sqrt(A)','i0','I_dot','Cus','Crc','Cis','Cic','Cuc','Crs','M0', 'OMEGA', 'e', 'Del_n','Toe','omega']
        dataset = pd.read_csv("5days.csv")
        dataset.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
        dataset = dataset.loc[dataset['PRN'] == int(sat_no)]
        df=dataset.copy()
        sss = df['SSS'][0]
        prn=df['PRN'][0]
        month = df['Epoch-Month'][0]
        year = df['Epoch-year'][0]
        sqrta=df['sqrt(A)'][0]
        i0=df['i0'][0]
        idot=df['I_dot'][0]

        dataset = dataset[(dataset['Epoch-Day'] >= fromhistory) & (dataset['Epoch-Day'] <= tohistory)]
        k=0
        m=[]
        for i in cols:
          col=i
          history = dataset[[col]].values
          predictions = list()
          try:
               for t in range(timestamp):
                  model = ARIMA(history, order=tuple(order[k]))
                  model_fit = model.fit(disp=0)
                  output = model_fit.forecast()
                  yhat = output[0]
                  predictions.append(yhat)
                  history = np.append(history, yhat)
                  if(k==0):
                      m.append([sss,prn,year,month,fromreq+int(t/20),(t%20)*2,0,0,sqrta,i0,idot])
                  m[t].extend(yhat)
          except (ValueError,np.linalg.LinAlgError):
              for t in range(timestamp):
                  model = ARIMA(history, order=(0, 0, 1))
                  model_fit = model.fit(disp=0)
                  output = model_fit.forecast()
                  yhat = output[0]
                  print(yhat)
                  predictions.append(yhat)
                  history = np.append(history, yhat)
                  if (k == 0):
                      m.append([sss,prn,year,month,fromreq+int(t/20),(t%20)*2,0,0,sqrta,i0,idot])
                  m[t].extend(yhat)
          k=k+1
          # plot
          plt.plot(predictions, color='red')
          plt.show()
          print(m)
          print('m0',len(m[0]))
          print('m1',len(m[1]))
        ga = pd.DataFrame(data=m, columns=head)
        ga.to_csv('hello.csv')
        return render_template("mod.html")

# Instantiate the solver
# ks = pyasl.MarkleyKESolver()

# Solves Kepler's Equation for a set
# of mean anomaly and eccentricity.
# Uses the algorithm presented by
# Markley 1995.
# M = 0.75
# e = 0.3
# print("Eccentric anomaly: ", ks.getE(M, e))




if __name__ == '__main__':
    app.run()

