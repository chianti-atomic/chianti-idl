
PRO WGFA_TIDY, FILEIN, FILEOUT, ENFILE=ENFILE, NO_INFO=NO_INFO, $
               wvl_decimal=wvl_decimal, quiet=quiet, $
               recompute_zero_wvl=recompute_zero_wvl, $
               text_output=text_output

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
;                    integer, this can be over-ridden. **This is now
;                    obsolete**. 
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
; OPTIONAL OUTPUTS:
;       Text_Output:  A string array containing the information
;                     that is usually printed directly to the IDL text
;                     window.
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
;       Ver.9, 14-Jan-2019, Peter Young
;                Now automatically sets the number of decimal places
;                for the wavelengths, so wvl_decimal is now obsolete;
;                added a comment to the output .wgfa stating when the
;                file was tidied; changed "close,lout" to
;                "free_lun,lout"; when applying the A-value cutoff I
;                do not take account of autoionizations or two-photon
;                transitions; introduced TEXT_OUTPUT= optional
;                output.
;       Ver.10, 6-Mar-2019, Peter Young
;                Changed how theoretical wavelengths are
;                calculated. If either of the levels does not have an
;                observed energy, then the wavelength is now
;                calculated only with the theoretical energies. 
;
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> wgfa_tidy, filein, fileout [, enfile=, /no_info, /quiet, text_output= ]'
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
 ; PRY, 6-Mar-2019
 ; I've modified this to reflect the recommendation given in
 ; CHIANTI Technical Report No. 10. That is, if either of the two
 ; energies is theoretical, then the wavelength must be derived from
 ; only the theoretical energies.
 ;
  k=where(en_l1 EQ -1. OR en_l2 EQ -1.,nk)
  IF nk NE 0 THEN BEGIN
    en_l1[k]=enth_l1[k]
    en_l2[k]=enth_l2[k]
    en_flag[k]=1b
  ENDIF 
  ;; k=where(en_l1 EQ -1.,nk)
  ;; IF nk NE 0 THEN BEGIN
  ;;   en_flag[k]=1b
  ;;   en_l1[k]=enth_l1[k]
  ;; ENDIF 
 ;
  ;; k=where(en_l2 EQ -1.,nk)
  ;; IF nk NE 0 THEN BEGIN
  ;;   en_flag[k]=1b
  ;;   en_l2[k]=enth_l2[k]
  ;; ENDIF
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

text_output=''

FOR i=1,n DO BEGIN
  k=where(l2 EQ i,nk)
  IF nk NE 0 THEN BEGIN
    avals=a[k]
    w=wvl[k]
   ;
   ; When applying the A-value cutoff I do not want to include
   ; autoionization or 2-photon rates.
   ;
    chck=where(w NE 0.,nchck)
    IF nchck GT 0 THEN BEGIN
      max_aval_i=max(avals[chck])
    ENDIF ELSE BEGIN
      max_aval_i=0.
    ENDELSE
   ;
    kk=where(avals GE max_aval_i*aval_cutoff OR w EQ 0.,nkk)
    kkx=where(avals LT max_aval_i*aval_cutoff AND w NE 0.,nkkx)
    IF nkkx NE 0 THEN BEGIN
      FOR j=0,nkkx-1 DO BEGIN
        ix=k[kkx[j]]
        text_output=[text_output,'   Transition ('+trim(l1[ix])+','+trim(l2[ix])+') -- removed']
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
          IF w EQ 0. THEN BEGIN
            print,'% WGFA_TIDY: Two levels have the same energy. Please check elvlc file.'
            print,'        Levels: ',l1[ind],l2[ind]
            return
          ENDIF 
          IF en_flag[ind] EQ 1 THEN w=-w
         ;
          transition_string=elvlc.data[l1[ind]-1].full_level+' - '+ $
                            elvlc.data[i-1].full_level
        ENDELSE
      ENDELSE
     ;
      IF keyword_set(no_info) THEN transition_string=''
     ;
     ; 14-Jan-2019, PRY
     ; Here I set different numbers of decimal places for the
     ; wavelengths based on the wavelength value.
     ;
      CASE 1 OF
        abs(w) LT 50.: format='(2i5,f15.4,2e15.3)'
        abs(w) GE 5e4 AND abs(w) LT 5e5: format='(2i5,f15.2,2e15.3)'
        abs(w) GE 5e5: format='(2i5,f15.2,2e15.3)'
        ELSE: format='(2i5,f15.3,2e15.3)'
      ENDCASE         
     ;
      datastr=string(format=format,l1[ind],l2[ind],w,gf[ind],a[ind])
      printf,lout,datastr+'   '+transition_string
    ENDFOR 
  ENDIF 
ENDFOR 

;
; Get a date-stamp
;
jd=systime(/julian,/utc)
mjd=jd-2400000.5d
mjd_str={ mjd: floor(mjd), time: (mjd-floor(mjd))*8.64d7 }
date_stamp=anytim2utc(/vms,mjd_str,/date)

PRINTF,lout,FORMAT='(i3)',-1
n=N_ELEMENTS(ref)
FOR i=0,n-1 DO PRINTF,lout,ref(i)
IF NOT keyword_set(no_comment) THEN BEGIN
  printf,lout,'%File processed with wgfa_tidy by '+getenv('USER')+' on '+date_stamp
ENDIF 
PRINTF,lout,FORMAT='(i3)',-1

free_lun,lout

;
; Create text summarizing results and put in text_output. Print text
; to screen if /quiet is not set.
;
IF ai_count NE 0 THEN BEGIN
  text_output=[text_output, $
               '% WGFA_TIDY:  there are '+trim(ai_count)+' auto-ionizing transitions in this file.']
ENDIF 
IF twop_count NE 0 THEN BEGIN
  text_output=[text_output, $
               '% WGFA_TIDY:  there are '+trim(twop_count)+' two-photon transitions in this file.']
ENDIF 
text_output=[text_output, $
             '% WGFA_TIDY:  number of transitions in original file: '+trim(ntrans), $
             '% WGFA_TIDY:  number of transitions in new file:      '+trim(tr_count) ]
;
text_output=text_output[1:*]
;
IF NOT keyword_set(quiet) THEN BEGIN
  n=n_elements(text_output)
  FOR i=0,n-1 DO print,text_output[i]
ENDIF


END
