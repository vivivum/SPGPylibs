pro read_cmilos_result,file,dimx,dimy,model

col = fltarr(12,dimx,dimy)
openr,1,file
readf,1,col
close,1

model = fltarr(dimx,dimy,12)
model[*,*,0] = col[5,*,*]
model[*,*,1] = col[2,*,*]
model[*,*,2] = col[8,*,*]
model[*,*,3] = col[6,*,*]  
model[*,*,4] = col[7,*,*]
model[*,*,5] = col[3,*,*]
model[*,*,6] = col[4,*,*]
model[*,*,7] = col[9,*,*]
model[*,*,8] = col[10,*,*]
model[*,*,9] = 0. ;macro 
model[*,*,10] = 0. ;alpha 
model[*,*,11] = col[11,*,*] ;chi2
delvar,col

; INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]

;                                 printf("%d\n",contador);
;                                 printf("%d\n",iter);
;                                 printf("%f\n",initModel.B);
;                                 printf("%f\n",initModel.gm);
;                                 printf("%f\n",initModel.az);
;                                 printf("%f \n",initModel.eta0);
;                                 printf("%f\n",initModel.dopp);
;                                 printf("%f\n",initModel.aa);
;                                 printf("%f\n",initModel.vlos); //km/s
;                                 //printf("alfa \t:%f\n",initModel.alfa); //stay light factor
;                                 printf("%f\n",initModel.S0);
;                                 printf("%f\n",initModel.S1);
;                                 printf("%.10e\n",chisqrf);

end