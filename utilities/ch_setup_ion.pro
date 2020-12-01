
FUNCTION ch_setup_ion, ions, wvlmin=wvlmin, wvlmax=wvlmax, ioneq_file=ioneq_file, $
                       noprot=noprot, path=path, radtemp=radtemp, rphot=rphot, $
                       abund_file=abund_file, quiet=quiet, noionrec=noionrec, $
                       obs_only=obs_only, index_wgfa=index_wgfa, no_auto=no_auto,$
                       no_rrec=no_rrec


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
;
;       QUIET:   If set, then information messages are not printed.
;       NOPROT:  If set, then proton rates are not read for the ion,
;                even if they exist.
;       NOIONREC: If set, then level-resolved ionization/recombination
;                rates are not read for the ion, even if they exist.
;       OBS_ONLY: If a wavelength check is performed (WVLMIN and/or
;                 WVLMAX), then the default is to consider observed
;                 and theoretical wavelengths. With /OBS_ONLY, only
;                 the observed wavelengths are considered.
;       NO_AUTO: If set, then the autoionization rates (contained in
;                the .auto file) are not read.
;       NO_RREC: If set, do not read the level-resolved radiative
;                recombination (RR) files.
;
; OUTPUTS:
;       A structure with the tags:
;       .gname   Name of the ion in CHIANTI format.
;       .jj      Array containing J-values for all levels.
;       .ecm     Array containing level energies. Energies are
;                observed if available otherwise they are
;                theoretical. 
;       .ecmth   Array containing theoretical energies.
;       .wvl     2D array containing wavelengths. No negative values.
;       .a_value 2D array containing A-values.
;       .wgfastr Structure containing output from read_wgfa_str.
;       .splstr  Structure containing output from read_scups.
;       .ioneq_file The name of the ion balance file.
;       .ip      Ionization potential in cm^-1 units.
;       .two_photon A structure with two tags containing the
;                two-photon rate and upper level. If data not
;                available, then it is a scalar with value -1.
;
;       The following tags will be included if the relevant data-sets
;       exist:
;       .prot_struc  Structure containing proton rate data.
;       .ionrec  Structure containing level-resolved ionization and
;                recombination data.
;       .abund_file  The name of an element abundance file.
;       .dilute  The radiation dilution factor (derived from RPHOT).
;       .radtemp The value of RADTEMP.
;       .autostr Structure containing autoionization rates from the
;                .auto files (same format as returned by read_auto). 
;       .rrec    for the level-resolved RR
; 
;       If a problem is found, then the integer -1 is returned.
;
; OPTIONAL OUTPUTS:
;       Index_Wgfa:  This is an index of the WGFA structure that picks
;                    out those transitions that satisfy the WVLMIN
;                    and/or WVLMAX conditions, and the A_VALUE NE 0
;                    condition. This index is used by the routine
;                    CH_SYNTHETIC routine (see "anylines" in this
;                    routine). 
;
; EXAMPLE:
;       IDL> output=ch_setup_ion('o_6')
;       IDL> output=ch_setup_ion('fe_13',wvlmin=170,wvlmax=210)
;       IDL> output=ch_setup_ion('si_12',path='/mydata/si_12')
;
; CALLS:
;       READ_WGFA_STR, CONVERTNAME, READ_ELVLC, READ_SCUPS,
;       READ_IONREC, R2W, READ_SPLUPS, READ_AUTO, CH_IP, read_rrlvl
;
; MODIFICATION HISTORY:
;      Ver.1, 28-Jun-2017, Peter Young
;         Modified from setup_ion.pro.
;      Ver.2, 9-Aug-2017, Peter Young
;         Added the elvlc structure to the output; implemented
;         /noionrec.
;      Ver.3, 25-Jan-2018, Peter Young
;         Added wgfastr to the output structure (needed this for
;         bb_rad_loss).
;      Ver.4, 31-Jan-2018, Peter Young
;         Fixed bug with wvlmax implementation; added /OBS_ONLY
;         keyword; removed call to read_wgfa2 (only call read_wgfa_str
;         now).
;      Ver.5, 4-Feb-2018, Peter Young
;         Added INDEX_WGFA optional output.
;      Ver.6, 3-Jun-2018, Peter Young
;         Reads the new .auto files (released with v.9) and puts them
;         in the structure.
;      Ver.7, 5-Jun-2018, Peter Young
;         Added ionization potential to output structure; added
;         information to header and comments.
;      Ver.8, 10 Oct 2018, Giulio Del Zanna (GDZ)
;         return -1 if no main data files (.wgfa, .elvlc, .scups) are
;         found.  
;      Ver.9, 12-Dec-2018, GDZ 
;         Added reading the new level-resolved RR files.
;      Ver.10, 12-Nov-2020, Peter Young
;         Added the tag 'two_photon' to the output. This is the
;         optional structure returned by read_wgfa_str.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output = ch_setup_ion ( ion_name [, wvlmin=, wvlmax=, ioneq_file=, '
  print,'                                   /noionrec, /noprot, /quiet, abund_file=, '
  print,'                                   rphot=, radtemp=, path=, index_wgfa= '
  print,'                                   /obs_only, /no_auto ])'
  return,-1
ENDIF

IF n_elements(ioneq_file) EQ 0 THEN BEGIN
  ioneq_file=!ioneq_file
ENDIF ELSE BEGIN
  chck=file_search(ioneq_file,count=count)
  IF count EQ 0 AND NOT keyword_set(quiet) THEN BEGIN
    print,'% CH_SETUP_ION: the specified IONEQ_FILE was not found. Returning...'
    return,-1
  ENDIF 
ENDELSE 

IF N_ELEMENTS(path) NE 0 THEN fname = concat_dir(path, ions) $
                         ELSE ion2filename,ions,fname
wname=fname+'.wgfa'
elvlcname=fname+'.elvlc'
upsname=fname+'.scups'
pname=fname+'.psplups'
autoname=fname+'.auto'

chck=file_search(wname,count=count1)
chck=file_search(elvlcname,count=count2)
chck=file_search(upsname,count=count3)


IF count1 eq 0 or count2 eq 0 or count3 eq 0 then begin 
if not keyword_set(quiet) then print, '% CH_SETUP_ION: no data files available for the ion ',ions
return, -1 
end 


convertname,ions,iz,ion

;
; Read .wgfa file.
; ---------------
read_wgfa_str,wname,wgfastr,two_photon=two_photon


;
; Use wvlmin and wvlmax to see if there are any lines in the requested
; wavelength range.
;
IF keyword_set(obs_only) THEN wvlchck=wgfastr.wvl ELSE wvlchck=abs(wgfastr.wvl)
;
CASE 1 OF
  n_elements(wvlmin) NE 0 AND n_elements(wvlmax) EQ 0: $
     k=where(wvlchck GE wvlmin AND wgfastr.aval NE 0.,nk)
  n_elements(wvlmin) EQ 0 AND n_elements(wvlmax) NE 0: $
     k=where(wvlchck LE wvlmax AND wgfastr.aval NE 0.,nk)
  n_elements(wvlmin) NE 0 AND n_elements(wvlmax) NE 0: $
     k=where(wvlchck LE wvlmax AND wvlchck GE wvlmin AND wgfastr.aval NE 0.,nk)
  ELSE: BEGIN
    nk=n_elements(wgfastr.wvl)
    k=lindgen(nk)
  ENDELSE 
ENDCASE
index_wgfa=k

IF nk EQ 0 THEN BEGIN
  index_wgfa=-1
  IF NOT keyword_set(quiet) THEN print,'% CH_SETUP_ION: no transitions found for '+trim(ions)+'. Returning...'
  return,-1
ENDIF 



;
; Process radiative data to get WVL and A_VALUE 2D arrays. Note that
; the routine checks for transitions with zero wavelengths (2-photon
; and autoionization rates) and adds these rates to the A-values.
; Also note that WVL is forced to be positive.
;
lvl1=wgfastr.lvl1
lvl2=wgfastr.lvl2
nlvls=max([lvl1,lvl2])
wvl=fltarr(nlvls,nlvls)
gf=fltarr(nlvls,nlvls)
a_value=fltarr(nlvls,nlvls)
wvl[lvl1-1,lvl2-1]= abs(wgfastr.wvl)
gf[lvl1-1,lvl2-1]= wgfastr.gf
a_value[lvl1-1,lvl2-1]=wgfastr.aval+wgfastr.auto


;
; Read .elvlc file to get JJ, ECM and ECMTH arrays
; ----------------
read_elvlc,elvlcname,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref,elvlcstr=elvlcstr
;
; The following sets the observed energy to be the theoretical energy
; in the case that the observed energy is missing (i.e., zero).
; 9-Nov-2015, PRY: I've added the line for eryd to be
; consistent, but eryd isn't used in the software.
;
; Recall that ecm=0 (not -1) for levels with no observed energy.
;
g=where(ecm EQ 0.)
IF max(g) GT 0 THEN BEGIN
  ecm[g]=ecmth[g]
  eryd[g]=erydth[g]
ENDIF 

;
; Read scups file to put collision data in SPLSTR
; ---------------
read_scups, upsname, splstr


;
; Get ionization potential in cm^-1 units.
; ---------------------------------------
ip=ch_ip(ions,/cm)

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
         wgfastr: wgfastr, $
         splstr:splstr, $
         ip: ip, $
         two_photon: two_photon, $
         ioneq_file: ioneq_file}


;
; Check for proton rates, and add structure to output
; ----------------------
chck=file_search(pname,count=count)
IF count NE 0 AND NOT keyword_set(noprot) THEN BEGIN
  IF NOT keyword_set(quiet) THEN BEGIN
    print,'% CH_SETUP_ION: proton rates added to output.'
  ENDIF 

  read_splups, pname, pstr, pref, /prot
  output=add_tag(output,pstr,'prot_struc')
ENDIF 


;
; Check for level-resolved ionization and recombination rates. 

 if file_exist(expand_path(fname+'.rrlvl')) and $ 
  file_exist(expand_path(fname+'.reclvl')) then begin 
 print, 'ERROR ! both  .reclvl and .rrlvl files found ! '
 return, -1
  end 

IF NOT keyword_set(no_rrec) THEN BEGIN
   rrec=read_rrlvl(fname, status)
  if status then begin 
     IF NOT keyword_set(quiet) THEN BEGIN
      print,'% CH_SETUP_ION: level-resolved radiative recombination rates added to output.'
    ENDIF 
    output=add_tag(output,rrec,'rrec')
 endif 
endif 


;
; Check for level-resolved ionization and recombination rates. If they
; exist, then create ionrec structure and add it as a tag to output.
;
IF NOT keyword_set(noionrec) THEN BEGIN
  read_ionrec,fname,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status
 ;
  IF status GT 0 THEN BEGIN
    IF NOT keyword_set(quiet) THEN BEGIN
      print,'% CH_SETUP_ION: level-resolved ionization & recombination rates added to output.'
    ENDIF 
   ;
    read_ioneq,ioneq_file,temp_all,ioneq,ioneq_ref
    IF ion GT 1 THEN ioneq_ionrec=reform(ioneq[*,iz-1,ion-2:ion])
    IF ion EQ 1 THEN ioneq_ionrec=reform(ioneq[*,iz-1,ion-1:ion])
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
    IF NOT keyword_set(quiet) THEN print,'% CH_SETUP_ION: the specified abundance file was not found, so it has not been added'
    print,'                                              to the output structure.'
  ENDELSE 
ENDIF 


;
; Add autoionization rates stored in the .auto files (if available).
;
IF file_exist(autoname) AND NOT keyword_set(no_auto) THEN BEGIN 
  IF NOT keyword_set(quiet) THEN BEGIN
    print,'% CH_SETUP_ION: autoionization rates added to output.'
  ENDIF 
 ;
  read_auto, autoname, autostr=autostr
  output=add_tag(output,autostr,'autostr')
ENDIF 



return,output

END
