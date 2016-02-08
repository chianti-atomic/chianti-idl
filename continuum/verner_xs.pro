
FUNCTION VERNER_XS, iz, ion, wvl, data=data

;+
; PROJECT:  CHIANTI
;
;      CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;      Astrophysical Plasmas. It is a collaborative project involving the Naval
;      Research Laboratory (USA), the University of Florence (Italy), the
;      University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME
;
;    VERNER_XS()
;
; EXPLANATION
;
;    Reads the Verner & Yakovlev (A&AS 109, 125, 1995) photoionization 
;    cross-section data and generates the values of the cross-section at 
;    the wavelengths WVL
;
; INPUTS
;
;    IZ     Atomic number of ion (e.g., 26 = Fe)
;
;    ION    Spectroscopic number of ion (e.g., 13 = XIII)
;
;    WVL    Wavelengths (in angstroms) for which cross-sections are 
;           required (1-D array).
;
; OPTIONAL INPUTS
;
;    DATA   By default VERNER_XS reads the Verner & Yakovlev data file when 
;           it is called. Through the DATA keyword the data array can be 
;           sent to VERNER_XS instead.
;
; OUTPUT
;
;    The photoionization cross-section for the ionization of the outer 
;    electron in units of mega-barns (Mb = 10^-18 cm^2) at the input 
;    wavelengths. E.g., for Fe XIII (ground configuration 
;    1s2.2s2.2p6.3s2.3p2) it is the cross-section for the removal of the 
;    3p electron.
;
; HISTORY
;
;    Ver.1, 24-Jul-2002, Peter Young
;-

IF n_elements(data) EQ 0 THEN BEGIN
  data=dblarr(10,465)
  dir=concat_dir(!xuvtop,'continuum')
  fname=concat_dir(dir,'verner_short.txt')
  openr,lun,fname,/get_lun
  readf,lun,data
  free_lun,lun
ENDIF

z=fix(data[0,*])
nel=fix(data[1,*])   ; no. of electrons

ions=z-nel+1

ind=where((iz EQ z) AND (ion EQ ions))

IF ind[0] EQ -1 THEN BEGIN
  print,'** No data for this ion! **'
  return, -1
ENDIF

IF n_elements(ind) GT 1 THEN BEGIN
  print,'** Something wrong with Verner data-file, please check **'
  return,-1
ENDIF

i=ind[0]

;
; convert wavelengths to energies in eV
;
wvl=double(wvl)
en=1d8/wvl/8065.54

m=n_elements(ind)

tot_sig=dblarr(n_elements(wvl))

n=data[2,i]
l=data[3,i]
e_th=data[4,i]
e_0=data[5,i]
sig_0=data[6,i]
y_a=data[7,i]
p=data[8,i]
y_w=data[9,i]
  
y=en/e_0

f_y= ((y-1)^2 +y_w^2) * $
     y^(-5.5-l+0.5*p) * $
     (1d0+sqrt(y/y_a))^(-p)

sig=f_y*sig_0

index=where(en LT e_th)
IF index[0] NE -1 THEN sig[index]=0d0

tot_sig=tot_sig+sig

return,tot_sig

END
