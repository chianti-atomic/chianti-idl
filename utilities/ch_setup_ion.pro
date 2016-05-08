
FUNCTION ch_setup_ion, ions, wvlmin=wvlmin, wvlmax=wvlmax, ioneq_file=ioneq_file, $
                       noprot=noprot, path=path, radtemp=radtemp, rphot=rphot, $
                       abund_file=abund_file, quiet=quiet, noionrec=noionrec


;+
; NAME:
;       CH_SETUP_ION
;
; PURPOSE:
;       A new version of setup_ion.pro that puts the atomic parameters
;       into a structure rather than common blocks. The structure can
;       be directly sent to the routine pop_solver for computing level
;       populations. 
;
; CATEGORY:
;       CHIANTI; data setup.
;
; CALLING SEQUENCE:
;	Result = CH_SETUP_ION( Ion_Name )
;
; INPUTS:
;	Ion_Name:  The name of an ion in CHIANTI format (e.g., 'o_6'
;	           for O VI).
;
; OPTIONAL INPUTS:
;	Wvlmin:  A wavelength in angstroms. If set, then the routine
;	         checks if the ion has any wavelengths above
;	         WVLMIN. If not, then the routine exits, and an empty
;	         output is returned.
;	Wvlmax:  A wavelength in angstroms. If set, then the routine
;	         checks if the ion has any wavelengths below
;	         WVLMAX. If not, then the routine exits, and an empty
;	         output is returned.
;       Ioneq_File: The name of an ionization equilibrium file, which
;                   is used to populate the level-resolved
;                   ionization/recombination data structure. If not
;                   specified then the file !ioneq_file is used.
;       Abund_File: The name of a CHIANTI format element abundance
;                   file. This would be used by pop_solver to compute
;                   the proton-to-electron ratio.
;       Radtemp: If photon excitation is included (by defining RPHOT),
;                then this input specifies the blackbody radiation
;                temperature in K. If not specified, then it is set to
;                6000 K.
;       Rphot:   Distance from the centre of the star in stellar radius units.
;                That is, RPHOT=1 corresponds to the star's
;                surface. If RPHOT is not specified, then photon
;                excitation will be switched off when pop_solver is
;                called.
;       Path:    This directly specifies the path where the
;                ion's data files are stored. If not set, then
;                the files are taken from the user's CHIANTI
;                distribution. 
;
; KEYWORD PARAMETERS:
;       QUIET:   If set, then information messages are not printed.
;       NOPROT:  If set, then proton rates are not read for the ion,
;                even if they exist.
;       NOIONREC: If set, then level-resolved ionization/recombination
;                rates are not read for the ion, even if they exist.
;
; OUTPUTS:
;       A structure with the tags:
;       .gname   Name of the ion in CHIANTI format.
;       .jj      Array containing J-values for all levels.
;       .ecm     Array containing observed energies.
;       .ecmth   Array containing theoretical energies.
;       .wvl     2D array containing wavelengths.
;       .a_value 2D array containing A-values.
;       .splstr  Structure containing output from read_scups.
;       .ioneq_file The name of the ion balance file.
;
;       The following tags will be included if the relevant data-sets
;       exist:
;       .prot_struc  Structure containing proton rate data.
;       .ionrec  Structure containing level-resolved ionization and
;                recombination data.
;       .abund_file  The name of an element abundance file.
;       .dilute  The radiation dilution factor (derived from RPHOT).
;       .radtemp The value of RADTEMP.
;
;       If a problem is found, then the integer -1 is returned.
;
; EXAMPLE:
;       IDL> output=ch_setup_ion('o_6')
;       IDL> output=ch_setup_ion('fe_13',wvlmin=170,wvlmax=210)
;       IDL> output=ch_setup_ion('si_12',path='/mydata/si_12')
;
; MODIFICATION HISTORY:
;      Ver.1, 28-Jun-2017, Peter Young
;         Modified from setup_ion.pro.
;      Ver.2, 9-Aug-2017, Peter Young
;         Added the elvlc structure to the output; implemented
;         /noionrec. 
;-



IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output = ch_setup_ion ( ion_name [, wvlmin=, wvlmax=, ioneq_file= '
  print,'                                   /noionrec, /noprot, /quiet, abund_file= '
  print,'                                   rphot=, radtemp=, path= ])'
  return,-1
ENDIF

IF n_elements(ioneq_file) EQ 0 THEN BEGIN
  ioneq_file=!ioneq_file
ENDIF ELSE BEGIN
  chck=file_search(ioneq_file,count=count)
  IF count EQ 0 AND NOT keyword_set(quiet) THEN BEGIN
    print,'%CH_SETUP_ION: the specified IONEQ_FILE was not found. Returning...'
    return,-1
  ENDIF 
ENDELSE 

IF N_ELEMENTS(path) NE 0 THEN fname = concat_dir(path, ions) $
                         ELSE ion2filename,ions,fname
wname=fname+'.wgfa'
elvlcname=fname+'.elvlc'
upsname=fname+'.scups'
pname=fname+'.psplups'

convertname,ions,iz,ion

;
; Read .wgfa file
;
read_wgfa2,wname,lvl1,lvl2,wvl1,gf1,a_value1,wgfaref

CASE 1 OF
  n_elements(wvlmin) NE 0 AND n_elements(wvlmax) EQ 0: $
     k=where(abs(wvl1) GE wvlmin AND a_value1 NE 0.,nk)
  n_elements(wvlmin) EQ 0 AND n_elements(wvlmax) NE 0: $
     k=where(abs(wvl1) LE wvlmax AND a_value1 NE 0.,nk)
  n_elements(wvlmin) NE 0 AND n_elements(wvlmax) EQ 0: $
     k=where(abs(wvl1) LE wvlmax AND abs(wvl1) GE wvlmin AND a_value1 NE 0.,nk)
  ELSE: nk=n_elements(lvl1)
ENDCASE 

IF nk EQ 0 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'%CH_SETUP_ION: no transitions found. Returning...'
  return,-1
ENDIF 



;
; Process radiative data to get WVL and A_VALUE 2D arrays. Note that
; the routine checks for transitions with zero wavelengths (2-photon
; and autoionization rates) and adds these rates to the A-values.
;
nlvls=max([lvl1,lvl2])
wvl=fltarr(nlvls,nlvls)
gf=fltarr(nlvls,nlvls)
a_value=fltarr(nlvls,nlvls)
;
ind1 = where(wvl1 EQ 0.)
ind2 = where(wvl1 NE 0.)
;
wvl[lvl1-1,lvl2-1]= abs(wvl1)
gf[lvl1-1,lvl2-1]=gf1
; 
IF ind1[0] NE -1 THEN BEGIN
  a_value[lvl1[ind1]-1,lvl2[ind1]-1] = a_value1[ind1]
  a_value[lvl1[ind2]-1,lvl2[ind2]-1] = $
     a_value[lvl1[ind2]-1,lvl2[ind2]-1] + a_value1[ind2]
ENDIF ELSE BEGIN
  a_value[lvl1[ind2]-1,lvl2[ind2]-1] = a_value1[ind2]
ENDELSE


;
; Read .elvlc file to get JJ, ECM and ECMTH arrays
;
read_elvlc,elvlcname,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref,elvlcstr=elvlcstr
;
; The following sets the observed energy to be the theoretical energy
; in the case that the observed energy is missing (i.e., zero).
; 9-Nov-2015, PRY: I've added the line for eryd to be
; consistent, but eryd isn't used in the software.
;
g=where(ecm EQ 0.)
IF max(g) GT 0 THEN BEGIN
  ecm[g]=ecmth[g]
  eryd[g]=erydth[g]
ENDIF 

;
; Read scups file to put collision data in SPLSTR
;
read_scups, upsname, splstr

;
; Initially create the "core" structure, i.e., containing the
; essential elements required by pop_solver.
; 
output={ gname: ions, $
         jj:jj,  $
         ecm:ecm, $
         ecmth:ecmth, $
         elvlcstr: elvlcstr, $
         wvl:wvl,  $
         a_value:a_value, $
         splstr:splstr, $
         ioneq_file: ioneq_file}


;
; Check for proton rates, and add structure to output
;
chck=file_search(pname,count=count)
IF count NE 0 AND NOT keyword_set(noprot) THEN BEGIN
  read_splups, pname, pstr, pref, /prot
  output=add_tag(output,pstr,'prot_struc')
ENDIF 



;
; Check for level-resolved ionization and recombination rates. If they
; exist, then create ionrec structure and add it as a tag to output.
;
IF NOT keyword_set(noionrec) THEN BEGIN 
  read_ionrec,fname,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status
 ;
  IF status GT 0 THEN BEGIN
    read_ioneq,ioneq_file,temp_all,ioneq,ioneq_ref
    IF ion GT 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-2:ion))
    IF ion eq 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-1:ion))
    ionrec = {rec:rec_rate, ci:ci_rate, temp:temp_ionrec, lev_up_rec:luprec, $
              lev_up_ci:lupci, status:status, ioneq:ioneq_ionrec, temp_ioneq:temp_all}
    output=add_tag(output,ionrec,'ionrec')
  ENDIF
ENDIF 

;
; Add the dilute and radtemp tags. Note that radtemp is ignored by
; pop_solver unless dilute is non-zero.
;
IF n_elements(rphot) NE 0 THEN BEGIN
  dilute=r2w(rphot)
  IF n_elements(radtemp) EQ 0 THEN radt=6d3 ELSE radt=double(radtemp)
  output=add_tag(output,dilute,'dilute')
  output=add_tag(output,radt,'radtemp')
ENDIF

;
; This adds the name of an abundance file. This will be used by
; pop_solver when computing the proton-to-electron ratio. If
; the tag doesn't exist, then pop_solver will just use the
; default abundance file. 
;
IF n_elements(abund_file) NE 0 THEN BEGIN
  chck=file_search(abund_file,count=count)
  IF count NE 0 THEN BEGIN
    output=add_tag(output,abund_file,'abund_file')
  ENDIF ELSE BEGIN
    IF NOT keyword_set(quiet) THEN print,'%CH_SETUP_ION: the specified abundance file was not found, so it has not been added'
    print,'                                              to the output structure.'
  ENDELSE 
ENDIF 


return,output

END
