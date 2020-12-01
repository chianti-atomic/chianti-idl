
PRO lorentz_test_3_top_lines

;+
; NAME:
;	LORENTZ_TEST_3_TOP_LINES
;
; PURPOSE:
;       This procedure computes the top 100 emission lines for
;       isothermal plasmas at three temperatures. This is the third of
;       the Lorentz tests for collisional ionization plasma codes. 
;
; CATEGORY:
;	Spectral modeling codes; output.
;
; CALLING SEQUENCE:
;       LORENTZ_TEST_3_TOP_LINES
;
; INPUTS:
;	None.
;
; OUTPUTS:
;       Creates text files containing the top 100 lines at the
;       specified temperatures.
;
; EXAMPLE:
;       IDL>  lorentz_test_3_top_lines
;
; MODIFICATION HISTORY:
;       Ver.1, 16-Aug-2016, Peter Young
;       Ver.2, 29-Jan-2019, Peter Young
;         Extended format string for the ion name as it was too
;         small.
;       Ver.3, 04-Nov-2020, Peter Young
;         Changed location of abundance file, to point to CHIANTI
;         abundance directory.
;-


wrange=[0.1,1000]   ; angstroms
density=1.0
volume=1e6    ; cm^3, equivalent to 1 m^3

temp=[1e6,6e6,4.642e7]
nt=n_elements(temp)

abund_dir=concat_dir(!xuvtop,'abundance')
abund_file=concat_dir(abund_dir,'proto_solar_2009_lodders.abund')
ioneq_file=!ioneq_file

read_abund,abund_file,ab,ref
nab=n_elements(ab)

FOR i=0,nt-1 DO BEGIN
 ;
 ; ch_synthetic takes as input the emission measure, which is
 ; N_e*N_H*V. For this test, V=10^6 cm^-3, and N_e=1.0. To compute N_H,
 ; we need to use the routine proton_dens
 ;
  n_h=proton_dens(alog10(temp[i]),/hydrogen,abund=abund_file,ioneq=ioneq_file)
  logem_iso=alog10(density^2*n_h*volume)

  ch_synthetic,wrange[0],wrange[1],output=output,density=density, $
               logt_iso=alog10(temp[i]), $
               logem_iso=logem_iso,ioneq=ioneq_file,/photon
 ;
 ; Multiply intensities by abundance
 ;
  FOR j=0,nab-1 DO BEGIN
    k=where(output.lines.iz EQ j+1,nk)
    IF nk GT 0 THEN output.lines[k].int=output.lines[k].int*ab[j]
  ENDFOR
 ;
 ; Get rid of steradian
 ;
  output.lines.int=output.lines.int*4.*!pi
 ;
 ; Now sort lines by intensity and extract top 100 lines
 ;
  k=reverse(sort(output.lines.int))
  lines=output.lines[k[0:99]]
 ;
 ; Open the output file
 ;
  tstr=trim(round(temp[i]/1e6))
  outfile='lorentz_test_3_chianti_'+tstr+'MK.txt'
  openw,lout,outfile,/get_lun
  printf,lout,'#Lorentz test: top 100 lines'
  printf,lout,'#Temperature: '+trim(string(format='(e12.3)',temp[i]))+' MK'
  printf,lout,'#Abundances: '+file_basename(abund_file)
  printf,lout,'#Density: '+trim(string(format='(e12.2)',density))+' cm^-3'
  printf,lout,'#Volume EM = N_e*N_H*V: '+trim(string(format='(e12.2)',10.^logem_iso))+' cm^-3'
  printf,lout,'#Column 1: line index [i3]'
  printf,lout,'#Column 2: wavelength (angstroms)  [f12.3]'
  printf,lout,'#Column 3: atomic number of emitting element  [i3]'
  printf,lout,'#Column 4: charge of emitting ion  [i3]'
  printf,lout,'#Column 5: log10 of intensity (photons s^-1)  [f10.3]'
  printf,lout,'#Column 6: ion name  [a11]'
  printf,lout,'#Column 7: transition information  [a40]'
  FOR j=0,99 DO BEGIN
    ion=strpad(lines[j].snote,11,fill=' ',/after)
    id=strpad(lines[j].ident,40,fill=' ',/after)
    printf,lout,format='(i3,f12.3,2i4,f10.3,3x,a11,a40)', $
           j+1,lines[j].wvl,lines[j].iz,lines[j].ion-1,alog10(lines[j].int),ion,id
  ENDFOR 
  free_lun,lout
  
ENDFOR 



END

