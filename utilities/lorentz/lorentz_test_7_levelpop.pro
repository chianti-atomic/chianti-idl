
PRO lorentz_test_7_levelpop, overwrite=overwrite, linelist=linelist, only_list=only_list


;+
; NAME:
;     LORENTZ_TEST_7_LEVELPOP
;
; PURPOSE:
;     Compute the results for Lorentz Test No. 7 "LevelPop", which
;     show the contributions of different atomic processes to key
;     X-ray lines. The list of lines is stored in CHIANTI and is
;     hardcoded to this routine.
;
; CATEGORY:
;     CHIANTI; Lorentz test.
;
; CALLING SEQUENCE:
;     LORENTZ_TEST_7_LEVELPOP
;
; INPUTS:
;     None.
;
; KEYWORD PARAMETERS:
;     OVERWRITE:  If the output file already exists, then this keyword
;                 forces it to be overwritten.
;     ONLY_LIST:  If set, then the routine only reads the line list
;                 (intended to be used in conjunction with
;                 LINELIST=). 
;
; OUTPUTS:
;     Creates a text file in the working directory with the results.
;
; OPTIONAL OUTPUTS:
;     LineList:  A structure containing the list of lines. The tags
;                are:
;                .ion  Ion name.
;                .lvl1  Lower level index.
;                .lvl2  Lower level index.
;                .wvl   Wavelength (angstroms).
;
; PROGRAMMING NOTES:
;     The quantities that are printed out are:
;      (i) Population out of level i by radiative process:
;           N_H*Ab*F*A_i*n_i
;          where N_H is number density of H, Ab is the element
;          abundance relative to H, F is the ionization fraction, A_i
;          is the rate out of level i to all levels, and n_i is the
;          population of level i relative to the ion.
;      (ii) Population out of level i by collisional process:
;           N_H*Ab*F*C_i*n_i*N_e
;           where C_i is the rate coefficient out of level i summed
;           over all levels, and N_e is the electron number density.
;
;     The final three columns in the file are dielectronic
;     recombination, inner shell ionization and inner shell
;     excitation, but these are just set to zero for CHIANTI. Cascades
;     following dielectronic capture are just part of the radiative
;     decay component. Inner shell processes (presumably to
;     autoionizing states that radiatively decay into the upper level)
;     are similarly not computed.
;
; MODIFICATION HISTORY:
;     Ver.1, 08-May-2019, Peter Young
;     Ver.2, 04-Nov-2020, Peter Young
;       Added /overwrite keyword. Changed the default location of the
;       output file, and the list of emission lines.
;     Ver.3, 05-Nov-2020, Peter Young
;       Fixed bug in previous version; added /only_list and linelist
;       inputs. 
;-


;
; The list of emission lines for which output is produced is stored
; in the lorentz directory of the CHIANTI IDL tree.
;
which,'lorentz_test_7_levelpop',/quiet,outfile=outfile
listdir=file_dirname(outfile)
listfile=concat_dir(listdir,'test7_lines.txt')


openr,lin,listfile,/get_lun
;
str={ion: '', lvl1: 0, lvl2: 0, wvl: 0.}
linelist=0
;
WHILE eof(lin) NE 1 DO BEGIN
  readf,lin,format='(a5,i4,i5,f9.0)',str
  IF n_tags(linelist) EQ 0 THEN linelist=str ELSE linelist=[linelist,str]
ENDWHILE 
;
free_lun,lin
;
IF keyword_set(only_list) THEN return
;
n=n_elements(linelist)


abundfile=concat_dir(!xuvtop,'abundance')
abundfile=concat_dir(abundfile,'proto_solar_2009_lodders.abund')
read_abund,abundfile,abund,ref


outfile='lorentz_test_7_chianti.txt'
chck=file_search(outfile,count=count)
IF count NE 0 AND NOT keyword_set(overwrite) THEN BEGIN
   print,'% LORENTZ_TEST_7_LEVELPOP: the output file already exists. Use /overwrite to overwrite it. Returning...'
   return
ENDIF 



temp=[1e6,6e6,4.642e7]
nt=n_elements(temp)

dens=[1.0,1e12]
nd=n_elements(dens)

openw,lout,outfile,/get_lun

printf,lout,'#Column 1: Index of emission line (see reference paper for list of lines and index).'
printf,lout,'#Column 2: Temperature (K) at which processes are calculated.'
printf,lout,'#Column 3: Electron number density (cm^-3) at which processes are calculated.'
printf,lout,'#Column 4: Sum of electron excitation rates (cm^-3 s^-1) into and out of level.'
printf,lout,'#Column 5: Sum of electron de-excitation rates (cm^-3 s^-1) into and out of level.'
printf,lout,'#Column 6: Sum of proton excitation rates (cm^-3 s^-1) into and out of level.'
printf,lout,'#Column 7: Sum of proton de-excitation rates (cm^-3 s^-1) into and out of level.'
printf,lout,'#Column 8: Sum of radiative decay rates (cm^-3 s^-1) into the level.'
printf,lout,'#Column 9: Sum of radiative decay rates (cm^-3 s^-1) out of the level.'
printf,lout,'#Column 10: Sum of radiative recombination rates (cm^-3 s^-1) into the level.'
printf,lout,'#Column 11: Sum of dielectronic recombination rates (cm^-3 s^-1) into the level.'
printf,lout,'#Column 12: Sum of inner shell ionization rates (cm^-3 s^-1) into the level.'
printf,lout,'#Column 13: Sum of inner shell excitation rates (cm^-3 s^-1) into the level.'
printf,lout,'#Derived using CHIANTI version '+ch_get_version()
printf,lout,'#Abundance file: '+file_basename(abundfile)
printf,lout,'#File created: '+systime()


;
; Note that collision rates from pop_processes are already multiplied
; by the no. density of the collider (electrons, protons).
;
FOR i=0,n-1 DO BEGIN
  FOR j=0,nt-1 DO BEGIN
    FOR k=0,nd-1 DO BEGIN
      convertname,trim(linelist[i].ion),iz,ion
      ab=abund[iz-1]
      ieq=get_ieq(temp[j],iz,ion,ioneq_name=!ioneq_file)
      nhne=proton_dens(alog10(temp[j]),abund=abundfile,ioneq=!ioneq_file,/hydrogen)
     ;
      pop_processes,trim(linelist[i].ion),level=linelist[i].lvl2, temp=temp[j], dens=dens[k], $
                    output=output,/quiet
     ;
      printf,lout,format='(i3,12e10.2)',i+1,temp[j],dens[k], $
             (output.in.e_exc+output.out.e_exc)*ab*ieq*nhne*dens[k], $
             (output.in.e_deexc+output.out.e_deexc)*ab*ieq*nhne*dens[k], $
             (output.in.p_exc+output.out.p_exc)*ab*ieq*nhne*dens[k], $
             (output.in.p_deexc+output.out.p_deexc)*ab*ieq*nhne*dens[k], $
             (output.in.rad_decay)*ab*ieq*nhne*dens[k], $
             (output.out.rad_decay)*ab*ieq*nhne*dens[k], $
             (output.in.rr)*ab*ieq*nhne*dens[k], $
             0., $    ; diel. recomb.
             0., $    ; inner shell ionization
             0.       ; inner shell excitation
             
    ENDFOR 
  ENDFOR 
ENDFOR

free_lun,lout

END
