
PRO WGFA_TIDY, FILEIN, FILEOUT, ENFILE=ENFILE, NO_INFO=NO_INFO, $
               wvl_decimal=wvl_decimal, quiet=quiet, $
               recompute_zero_wvl=recompute_zero_wvl

;+
; NAME:
;       WGFA_TIDY
;
; PURPOSE:
;	To reduce the size of Chianti .wgfa files
;
; EXPLANATION:
;       This routine takes an existing CHIANTI .wgfa file and
;       re-writes it, removing any weak A-values. Optionally the
;       routine will take the energies from the CHIANTI .elvlc file
;       and re-compute the wavelengths.
;
;       Weak A-values are defined as those that are more than a
;       factor 10^7 weaker than the strongest transition coming from
;       that level.
;
;       The wavelengths are computed from the elvlc file using the
;       observed energies. If one or both of the levels does not have
;       an observed energy then theoretical energy for the level(s) is
;       used.
;
;       The routine checks for transitions with a zero wavelength in
;       the original .wgfa file and assumes that these are two-photon
;       or autoionizing transitions. The wavelengths are not
;       re-computed. 
;
; INPUTS:
;	FILEIN	The original Chianti .wgfa file
;	FILEOUT	The name of the revised Chianti file
;
; OPTIONAL INPUTS:
;       ENFILE:  Allows the input of a specified .elvlc file.
;
;       WVL_DECIMAL: By default the routine prints wavelengths to 3
;                    decimal places. Specifying wvl_decimal to be an
;                    integer, this can be over-ridden.
;
; KEYWORD PARAMETERS:
;       NO_INFO: By default the routine appends transition information
;                to the end of each data line. Setting /NO_INFO
;                suppresses this.
;       QUIET:   If set, then no information messages will be printed.
;
;       RECOMPUTE_ZERO_WVL: By default the routine will not recompute
;                           the wavelengths of transitions that have a
;                           zero wavelength (since these are assumed
;                           to be two-photons or
;                           autoionizations). Setting this keyword
;                           forces the wavelengths to be re-computed.
;
; CALLS
;       CHIANTI:  READ_WGFA2, READ_ELVLC, CH_IP
;
; EXAMPLE
;
;	WGFA_TIDY, !xuvtop+'/fe/fe_14/fe_14.wgfa', 'fe_14.wgfa_new'
;
; HISTORY
;
;	Ver.1, 9-Feb-99, PRY
;       Ver.2, 24-Aug-00, PRY
;       Ver.3, 11-Dec-00, PRY
;                Added /USE_RYD keyword.
;       Ver.4, 5-Nov-02, PRY
;                Disabled /USE_RYD option following change to 
;                read_best_energs.pro 
;       Ver.5, 28-Jan-2008, PRY
;                Changed A-value cutoff from 1e-5 to 1e-7.
;       Ver.6, 12-Feb-2009, PRY
;                Added /get_lun in call to openw, and removed direct
;                lun specification.
;       Ver.7, 5-Sep-2012, PRY
;                Added WVL_DECIMAL= keyword; updated header.
;       Ver.8, 14-Jul-2016, Peter Young
;                Routine re-written to handle the new elvlc files
;                distributed with CHIANTI 8.
;
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> wgfa_tidy, filein, fileout [, enfile=, /no_info, /quiet '
  print,'                          wvl_decimal= ]'
  return
ENDIF

;
; Print format for output file
;
IF n_elements(wvl_decimal) EQ 0 THEN wvl_decimal=3
format='(2i5,f15.'+trim(wvl_decimal)+',2e15.3)'



;
; read_wgfa2 puts the data into 1D arrays.
;
read_wgfa2,filein,l1,l2,wvl,gf,a,ref

ntrans=n_elements(l1)

;
; Read the elvlc file if it has been input.
;
IF N_ELEMENTS(enfile) NE 0 THEN BEGIN
  read_elvlc,enfile,lev,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth, $
             refen, elvlc=elvlc

  en_l1=elvlc.data[l1-1].obs_energy
  en_l2=elvlc.data[l2-1].obs_energy
  enth_l1=elvlc.data[l1-1].theory_energy
  enth_l2=elvlc.data[l2-1].theory_energy
  en_flag=bytarr(ntrans)
 ;
  k=where(en_l1 EQ -1.,nk)
  IF nk NE 0 THEN BEGIN
    en_flag[k]=1b
    en_l1[k]=enth_l1[k]
  ENDIF 
 ;
  k=where(en_l2 EQ -1.,nk)
  IF nk NE 0 THEN BEGIN
    en_flag[k]=1b
    en_l2[k]=enth_l2[k]
  ENDIF
 ;
 ;
 ; Get ionization potential in order to check for autoionizing levels
 ;
  ip=ch_ip(elvlc.info.ion_name)
ENDIF



;
; If max_aval is the maximum A-value from a given level, then only
; transitions from this level that have A-values >
; max_aval*aval_cutoff will be included in the output file.
;
aval_cutoff=1d-7


openw,lout,fileout,/get_lun

n=max([l1,l2])

ai_count=0
twop_count=0
tr_count=0

FOR i=1,n DO BEGIN
  k=where(l2 EQ i,nk)
  IF nk NE 0 THEN BEGIN
    avals=a[k]
    max_aval_i=max(avals)
    kk=where(avals GE max_aval_i*aval_cutoff,nkk)
    kkx=where(avals LT max_aval_i*aval_cutoff,nkkx)
    IF nkkx NE 0 AND NOT keyword_set(quiet) THEN BEGIN
      FOR j=0,nkkx-1 DO BEGIN
        ix=k[kkx[j]]
        print,'   Transition ('+trim(l1[ix])+','+trim(l2[ix])+') -- removed'
      ENDFOR 
    ENDIF 
    tr_count=tr_count+nkk
    FOR j=0,nkk-1 DO BEGIN
      ind=k[kk[j]]
      w=wvl[ind]
     ;
     ; Check if wavelength in existing file is zero - this implies the
     ; transition is a two-photon or autoionizing rate.
     ;
      IF w EQ 0. AND NOT keyword_set(recompute_zero_wvl) THEN BEGIN
        IF en_l2[ind] GT ip THEN BEGIN 
          transition_string=elvlc.data[i-1].full_level+' [autoionizing rate]'
          ai_count=ai_count+1
        ENDIF ELSE BEGIN
          transition_string=elvlc.data[l1[ind]-1].full_level+' - '+ $
                            elvlc.data[i-1].full_level+' [two-photon]'
          twop_count=twop_count+1
        ENDELSE 
      ENDIF ELSE BEGIN 
        IF n_elements(enfile) EQ 0 THEN BEGIN
          transition_string=''
        ENDIF ELSE BEGIN 
          w=1d8/abs(en_l2[ind]-en_l1[ind])
          IF en_flag[ind] EQ 1 THEN w=-w
         ;
          transition_string=elvlc.data[l1[ind]-1].full_level+' - '+ $
                            elvlc.data[i-1].full_level
        ENDELSE
      ENDELSE
     ;
      IF keyword_set(no_info) THEN transition_string=''
     ;
      datastr=string(format=format,l1[ind],l2[ind],w,gf[ind],a[ind])
      printf,lout,datastr+'   '+transition_string
    ENDFOR 
  ENDIF 
ENDFOR 


PRINTF,lout,FORMAT='(i3)',-1
n=N_ELEMENTS(ref)
FOR i=0,n-1 DO PRINTF,lout,ref(i)
PRINTF,lout,FORMAT='(i3)',-1

CLOSE,lout


IF ai_count NE 0 AND NOT keyword_set(quiet) THEN BEGIN
  print,'% WGFA_TIDY:  there are '+trim(ai_count)+' auto-ionizing transitions in this file.'
ENDIF 
IF twop_count NE 0 AND NOT keyword_set(quiet) THEN BEGIN
  print,'% WGFA_TIDY:  there are '+trim(twop_count)+' two-photon transitions in this file.'
ENDIF 
IF NOT keyword_set(quiet) THEN BEGIN
  print,'% WGFA_TIDY:  number of transitions in original file: '+trim(ntrans)
  print,'% WGFA_TIDY:  number of transitions in new file:      '+trim(tr_count)
ENDIF

END
