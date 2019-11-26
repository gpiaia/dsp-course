#!/usr/bin/env python3.4
# -*- coding: utf-8 -*-

"""
Created on Wed Jul 17 19:07:37 2019

@author: Cesar Crovato
"""

import numpy as np
from scipy.fftpack import fft
from scipy.signal import hanning

class fft_interpol:
    def __init__(self):
         print('Init Interpolation')

    def calculo_espectro(fs,maximo_numeroComponentes,limiar,janela , databuf):

        def detecta_picos(maximo_numeroPicos,limiar,Modulos):


            index = 0
            MaiorP=0
            MaiorIndex=0

            n1 = Modulos[0]
            n2 = Modulos[1]
            n3 = Modulos[2]
            n4 = Modulos[3]
            Saida=[]

            for k in range(3,len(Modulos)):
                n4=Modulos[k]
                if ((n2>n1) and (n2>n3) and (n2>limiar)):
                    Saida.append(k-2)
                    index=index+1
                    if (Modulos[Saida[index-1]]>MaiorP):
                        MaiorIndex=Saida[index-1]
                        MaiorP=Modulos[MaiorIndex]




                if ((n2>n1) and (n2==n3) and (n3>n4) and (n2>limiar)):
                    Saida.append(k-2)
                    index=index+1
                    if (Modulos[Saida[index-1]]>MaiorP):
                        MaiorIndex=Saida[index-1]
                        MaiorP=Modulos[MaiorIndex]
                n1=n2
                n2=n3
                n3=n4
                if (len(Saida)>=maximo_numeroPicos):
                    break
            return MaiorP,MaiorIndex,Saida



        def signum(val):
            if abs(val)<1e-10:
                return 0
            else:
                if (val > 0):
                    return  1
                else:
                    if (val < 0):
                        return -1


        def Interpola(fs,maximo_numeroComponentes,limiar,janela,Modulos,Angulos):

            PicoMax,IndPicoMax,indpicos=detecta_picos(maximo_numeroComponentes,limiar,Modulos)
            Mod_interp=[]
            Ang_interp=[]
            Freq_interp=[]

            if PicoMax==0:
                return 0,0,0,[0],[0],[0]
            else:
                num_pts_espectrais=len(Modulos)

                for i in range (0,len(indpicos)):
                    index=indpicos[i]
                    if janela=='retangular':
                        print("os resultados são piores que hanning não, usar retangular")
                        mod_ant =  (Modulos[index-1])
                        mod_pos =  (Modulos[index+1])
                        mod =     (Modulos[index])
                        delta1=(mod_pos+mod_ant)/(2.0*mod-mod_ant+mod_pos)
                        delta2=(mod_pos)/(mod+mod_pos)
                        delta = 0.5*(delta1+delta2)
                        abs_delta = np.abs(delta)
                        if(mod_pos>=mod_ant):
                            delta = -abs_delta
                        else:
                            delta = abs_delta

                        f=index-delta
                        f=f*fs/(num_pts_espectrais*2)

                        if (abs_delta<=0.0000001):
                            Amp=mod
                            Ang=Angulos[index]
                        else:
                            Amp=np.pi*mod*delta/np.sin(np.pi*delta)
                            #//Ang= (Angulos[index])+delta*np.pi*(num_pts_espectrais*2-1)/(num_pts_espectrais*2)  //em radianos
                            Ang=Angulos[index]


                    if janela=='hanning':

                        try:
                            sig_delta=signum(Modulos[index+1]-Modulos[index-1])
                            mod =  2*Modulos[index]
                            alfa=Modulos[index]/Modulos[index+sig_delta]
                            mod_delta=(2-alfa)/(1+alfa)
                            delta=mod_delta*sig_delta


                        except:
                            delta=0.0

                        f=index+delta
                        f=f*fs/(num_pts_espectrais*2)
                        #stdoutmutex.acquire()
                        #print("freq %f",f)
                        #stdoutmutex.release()

                        if (delta==0.0):
                             Amp= mod
                             Ang=Angulos[index]
                        else:
                            try:
                                Amp=np.pi*mod*delta*(1-(delta*delta))/np.sin(np.pi*delta)
                            except:
                                Amp= mod
                            Ang=Angulos[index]-delta*np.pi
                            #Ang=Angulos[index]


                 #   if janela=='hanning':
                 #       Amp=Amp/2

                    if abs(Ang)<1e-5:
                        Ang=0
                    Mod_interp.append(Amp)
                    Ang_interp.append(Ang)
                    Freq_interp.append(f)

                    if index==IndPicoMax:
                        Mod_f0=Amp
                        Ang_f0=Ang
                        Freq_f0=f


                return Mod_f0,Ang_f0,Freq_f0,Mod_interp,Ang_interp,Freq_interp

        def rms(x,n)   :
            #se x é espectro, então n=2     se x é dominio tempo, então n deve ser len(x)
            return np.sqrt(np.sum(np.array(x)*np.array(x))/n)


        def cosseno(Amp,Fase,Offset,Freq,Np,Fs):
            incremento=1/Fs
            ti=0
            tf=incremento*(Np-1)
            t = np.linspace(ti, tf, num=Np)
            y = Offset + Amp*np.cos(2*np.pi*Freq*t+Fase*np.pi/180)
            return y

        N=len(databuf)
        if janela=='hanning':
            databuf_comjanelamento=databuf* hanning(len(databuf),sym=False)

        if janela=='retangular':
            databuf_comjanelamento=databuf

        yf = fft(databuf_comjanelamento)
        yf=yf[0:int(N/2)]

        Modulos=2.0/N * np.abs(yf)
        Angulos=np.angle(yf, deg=False)
        Freq = np.linspace(0, (N/2)*fs/N,num=N/2)


        Mod_f0,Ang_f0,Freq_f0,Mod_interp,Ang_interp,Freq_interp= Interpola(fs,maximo_numeroComponentes,limiar,janela,Modulos,Angulos)
        Modulos[np.where(np.array(Modulos)<=limiar)]=0 #zeramos de acordo com o limiar

        return Mod_f0,Ang_f0,Freq_f0,Mod_interp,Ang_interp,Freq_interp,Modulos,Angulos,Freq
