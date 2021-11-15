pro tachometer, landa_ref, deltaV, iprofile, vlos

;assuming first wavelength is the continuum point
; deltaV = sampling in volts
; vlos is given in km/s Volts/A
; Using the tunning constant it goes into km / s 

index = 1

Vlight = 2.99792458d18 ;A/s

tan_term = atan(iprofile(index)+iprofile(index+1)-iprofile(index+3)-iprofile(index+4),iprofile(index)-iprofile(index+1)-iprofile(index+3)+iprofile(index+4) )

vlos = 2. * Vlight * deltaV / !dpi / landa_ref * tan_term * 1.d-13

END

pro test

a=readfits('solo_L0_phi-fdt-ilam_0669785848_V202104090820C_0143230503.fits')

vl = fltarr(768, 768)

for i=0,768-1 do begin
    for j=0,768-1 do begin 
        tachometer, 6173.0 , 100., reform(a[i,j,*]), vlos
        vl[i,j] = vlos 
    endfor
endfor
stop
end

pro test2

DEVICE, GET_DECOMPOSED=old_decomposed
DEVICE, DECOMPOSED=0
WINDOW, 1, XSIZE=1200, YSIZE=1200
userlct,coltab=7,/full

f=file_search('helio/*fits')
vl = fltarr(768, 768, n_elements(f))
for k=0,n_elements(f) - 1 do begin

    print,k,n_elements(f)-1
    a=readfits(f[k])
    for i=0,768-1 do begin
        for j=0,768-1 do begin 
            tachometer, 6172.7 , 200., reform(a[i,j,*]), vlos
            vl[i,j,k] = vlos 
        endfor
    endfor
    offset = median(vl[768/2-100:768/2+100,768/2-100:768/2+100,k])
    tvframe3,(vl[*,*,k]- offset)<6000>(-6000),title='Velocity [km/s Volts/A]',chars=2,charthick=2

    WRITE_PNG, 'dummy/image'+strtrim(string(k),1)+'.png', TVRD(/TRUE),XRESOLUTION=300,YRESOLUTION=300
endfor

end