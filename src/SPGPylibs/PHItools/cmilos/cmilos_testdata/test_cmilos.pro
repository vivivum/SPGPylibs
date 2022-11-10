; data comes from /Users/orozco/Dropbox (IdAdA)/TRABAJO/Codigo_Compresion/Compresion_DHernanez/datos_crisp
; fits_CRISP_1_6wl_binning3.txt
; CRISP_Sol_calma_001_EXP.txt

binning = 3.
tamX = nint(883./binning);
tamY = nint(894./binning);
dimmXY = float(tamx)*float(tamy)

IIi=fltarr(dimmXY,6)
QQq=fltarr(dimmXY,6)
UUu=fltarr(dimmXY,6)
VVv=fltarr(dimmXY,6)
wave_axis=fltarr(6)
rr = fltarr(5)

openr,1,'test_data.txt'&$
for i=0,dimmXY-1 do begin  &$
   for j=0,5 do begin  &$
      readf,1,rr  &$
      iii[i,j] = rr[1]  &$
      qqq[i,j] = rr[2]  &$
      uuu[i,j] = rr[3]  &$
      vvv[i,j] = rr[4]  &$
      wave_axis[j] = rr[0]  &$
   endfor  &$
   endfor
close,1

iii = reform(iii,tamX,tamY,6.)
qqq = reform(qqq,tamX,tamY,6.)
uuu = reform(uuu,tamX,tamY,6.)
vvv = reform(vvv,tamX,tamY,6.)

data=[[[iii]],[[qqq]],[[uuu]],[[vvv]]]
data = reform(data,tamX, tamY, 6,4)

;setpscf,filename='SPOT_CRISP.ps'
!p.multi=[0,6,4]
for i=0,5 do tvframe3,data[*,*,i,0],charsize=0.6,/bar,title='Wavelength '+text(i+1),ytitle='Stokes I'
for i=0,5 do tvframe3,data[*,*,i,1]<0.05>(-0.05),charsize=0.6,/bar,ytitle='Stokes Q'
for i=0,5 do tvframe3,data[*,*,i,2]<0.05>(-0.05),charsize=0.6,/bar,ytitle='Stokes U'
for i=0,5 do tvframe3,data[*,*,i,3]<0.1>(-0.1),charsize=0.6,/bar,ytitle='Stokes V'
;endpscf

save,filename='test_data.save',wave_axis,data

;CMILOS INVERSION

spawn,'.././milos 6 15 0 0 test_data.txt > output.txt'

;Leo resultados cmilos y los pongo en formato IDL_milos

;FIXME COMPARO LOS RESULTADOS
read_cmilos_result,'output.txt',294,298,cmodel


;***********
;IDL VERSION
;***********

init_milos,'6173',wl
smodel=fltarr(tamx,tamy,11)
imodel = [3,400.,0.01,0.025,1.,30.,120.,0.15,0.85,0.,0.]
for x = 0,tamx-1 do begin    &$
  for y = 0,tamy-1 do begin       &$
    perfil = reform(data[x,y,*,*])     &$
    model = imodel   &$
    milos, Wl, lll, model, perfil, chisqr=chisqr,yfit=yfit,$
      fix=[1,1,1,1,1,1,1,1,1,0,0],/inversion,miter=15,weight = [1,10,10,4],$
      ilambda=10,toplim=1D-12,numerical= 0,iter_info = iter_info,/quiet    &$
    smodel[x,y,*] = model   &$
  endfor   &$
  print,x,tamx &$
endfor

plot,smodel(*,*,2),cmodel(*,*,2),psym=3,xrange=[-3,3],yrange=[-3,3]

save,filename='test.sav',smodel













       ;inversiones

       B=fltarr(dimmXY)
       G=fltarr(dimmXY)
       A=fltarr(dimmXY)
       V=fltarr(dimmXY)
       rr = fltarr(10)

       openr,1,'crisp_milos_IDL_6wl_binning3.txt'&$
       for i=0,dimmXY-1 do begin  &$
             readf,1,rr  &$
             B[i] = rr[1]  &$
             G[i] = rr[5]  &$
             A[i] = rr[6]  &$
             V[i] = rr[2]  &$
          endfor
       close,1


       b = reform(b,tamY,tamX)
       g = reform(g,tamY,tamX)
       a = reform(a,tamY,tamX)
       v = reform(v,tamy,tamX)

       setpscf,filename='SPOT_CRISP_PARAMS.ps'
       !p.multi=[0,4,1]
       loadct,4
       tvframe3,b[*,*]<2500>0,charsize=0.6,title='Field strength [G] '
       tvframe3,g[*,*]<180>(0),charsize=0.6,title='Field inclination [degree] '
       tvframe3,a[*,*]<180>(0),charsize=0.6,title='Field azimuth [degree] '
       ;loadct,70
       ;loadcti
       ;!p.background=254
        tvframe3,v[*,*]<2>(-2),charsize=0.6,title='Doppler velocity [km/s] '
       endpscf
loadct,0


       ;QS

       binning = 1.
       tamX = nint(901./binning);
       tamY = nint(901./binning);
       dimmXY = float(tamx)*float(tamy)

              IIi=fltarr(dimmXY,6)
              QQq=fltarr(dimmXY,6)
              UUu=fltarr(dimmXY,6)
              VVv=fltarr(dimmXY,6)
              LLl=fltarr(6)
              rr = fltarr(5)

              openr,1,'CRISP_Sol_calma_001_EXP.txt'&$
              for i=0,dimmXY-1 do begin  &$
                 for j=0,5 do begin  &$
                    readf,1,rr  &$
                    iii[i,j] = rr[1]  &$
                    qqq[i,j] = rr[2]  &$
                    uuu[i,j] = rr[3]  &$
                    vvv[i,j] = rr[4]  &$
                    lll[j] = rr[0]  &$
                 endfor  &$
                 endfor
              close,1

              iii = reform(iii,tamX,tamY,6.)
              qqq = reform(qqq,tamX,tamY,6.)
              uuu = reform(uuu,tamX,tamY,6.)
              vvv = reform(vvv,tamX,tamY,6.)

              data=[[[iii]],[[qqq]],[[uuu]],[[vvv]]]
              data = reform(data,tamX, tamY, 6,4)

              setpscf,filename='QS_CRISP.ps'
              !p.multi=[0,6,4]
              for i=0,5 do tvframe3,data[100:800,100:800,i,0],charsize=0.6,/bar,title='Wavelength '+text(i+1),ytitle='Stokes I'
              for i=0,5 do tvframe3,data[100:800,100:800,i,1]<0.005>(-0.005),charsize=0.6,/bar,ytitle='Stokes Q'
              for i=0,5 do tvframe3,data[100:800,100:800,i,2]<0.005>(-0.005),charsize=0.6,/bar,ytitle='Stokes U'
              for i=0,5 do tvframe3,data[100:800,100:800,i,3]<0.01>(-0.01),charsize=0.6,/bar,ytitle='Stokes V'
              endpscf

;       crisp_milos_IDL_6wl_binning3.txt

tamx = 701.
tamy = 701.
tamxy = float(tamx)*float(tamy)
B=fltarr(tamxy)
G=fltarr(tamxy)
A=fltarr(tamxy)
V=fltarr(tamxy)
rr = fltarr(10)

openr,1,'crisp_qs_milos_idl_dos.txt'  &$
for i=0,tamxy-1 do begin &$
  readf,1,rr  &$
  B[i] = rr[1]  &$
  G[i] = rr[5]  &$
  A[i] = rr[6]  &$
  V[i] = rr[2]  &$
endfor
close,1

b = reform(b,tamY,tamX)
g = reform(g,tamY,tamX)
a = reform(a,tamY,tamX)
v = reform(v,tamy,tamX)

setpscf,filename='QS_CRISP_PARAMS.ps'
!p.multi=[0,4,1]
loadct,4
tvframe3,b[*,*]<2500>0,charsize=0.6,title='Field strength [G] '
tvframe3,g[*,*]<180>(0),charsize=0.6,title='Field inclination [degree] '
tvframe3,a[*,*]<180>(0),charsize=0.6,title='Field azimuth [degree] '
;loadct,70
;loadcti
;!p.background=254
 tvframe3,v[*,*]<2>(-2),charsize=0.6,title='Doppler velocity [km/s] '
endpscf
loadct,0




;DATOS BRUTOS SOL EN CALMA   CRISP_Sol_calma_001_EXP.txt
; 5 columnas 1 longitud de onda 2,3,4,5 IQUV
; 6 longitudes de onda

;PARAMETROS FISICOS


;crisp_qs_milos_idl_dos.txt
;PERO SOLO ESTE TROZO (100:800,100:800)
;10 parametros: 2 -> campo, 3 -> velocidad, 6 -> inclinacion, 7-> azimut


;LEO DATOS DE HMI AHORA


files = file_search('../hmi_images/hmi.ME_720s_fd10_2015.08.24_073000_TAI-2015.08.24_103000_TAI/*fits')

read_sdo,files(0),a,azi
read_sdo,files(1),a,fld
read_sdo,files(2),a,inc
read_sdo,files(3),a,vlo

!p.background=254/2
             setpscf,filename='HMI.ps'
       !p.multi=[0,4,1]
       loadct,4
		tvlct,rr,gg,bb,/get
		     rr(0) = 254
		     gg(0) = 254
		     bb(0) = 254
		tvlct,rr,gg,bb
       tvframe3,fld[*,*]<1000>0,charsize=0.6,title='Field strength [G] ',txcolor=1;,charcolor==250
       tvframe3,inc[*,*]<180>(0),charsize=0.6,title='Field inclination [degree] ',txcolor=1;,charcolor==250
       tvframe3,azi[*,*]<180>(0),charsize=0.6,title='Field azimuth [degree] ',txcolor=1;,charcolor=250
       ;loadct,70
       ;loadcti
       ;!p.background=254
        tvframe3,vlo[*,*]/100000.<2.5>(-2.5),charsize=0.6,title='Doppler velocity [km/s] ',txcolor=5
              endpscf



end
