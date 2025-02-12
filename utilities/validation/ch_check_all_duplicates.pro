

PRO ch_check_all_duplicates, outfile=outfile, quiet=quiet

;+
; NAME:
;     CH_CHECK_ALL_DUPLICATES
;
; PURPOSE:
;     Goes through the complete set of CHIANTI ions and checks for
;     duplicated transitions in the wgfa, scups and rrlvl files.
;
; CATEGORY:
;     CHIANTI; validity.
;
; CALLING SEQUENCE:
;     CH_CHECK_ALL_DUPLICATES
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     Outfile:  The name of a file to send the list of duplicated
;               transitions. If not set, they are send to
;               ch_duplicates.txt in the working directory.
;	
; KEYWORD PARAMETERS:
;     QUIET:  If set, then no information is written to the IDL window.
;
; OUTPUTS:
;     The list of problem files is printed to the IDL window (unless /quiet
;     set), and the list of problem transitions is written to
;     ch_duplicates.txt (or OUTFILE).
;
; CALLS:
;     CH_READ_LIST_IONS, CH_DUPLICATE_TRANSITIONS, CONVERTNAME,
;     ZION2FILENAME, Z2ELEMENT
;
; EXAMPLE:
;     IDL> ch_check_all_duplicates
;
; MODIFICATION HISTORY:
;     Ver.1, 07-Feb-2025, Peter Young
;-


mlist=ch_read_list_ions(count=n)
ions=mlist.list_ions

IF n_elements(outfile) EQ 0 THEN outfile='ch_duplicates.txt'

openw,lout,outfile,/get_lun

problem_files=''
problem_ions=''

file_ext=['.wgfa','.scups','.rrlvl']
n_ext=n_elements(file_ext)

progress,0.0,/reset,label='Progress (%)'
FOR i=0,n-1 DO BEGIN
  convertname,ions[i],iz,ion,diel=diel
  zion2filename,iz,ion,fname,diel=diel
  FOR j=0,n_ext-1 DO BEGIN 
    d=ch_duplicate_transitions(fname+file_ext[j],status=status,/quiet)
    IF status EQ 1 THEN BEGIN
      problem_files=[problem_files,file_basename(fname+file_ext[j])]
      problem_ions=[problem_ions,ions[i]]
      nt=n_elements(d)
      printf,lout,fname+file_ext[j]
      FOR k=0,nt-1 DO BEGIN
        trans=trim(d[k].lvl1)+'-'+trim(d[k].lvl2)
        trans=strpad(trans,10,/after,fill=' ')
        printf,lout,format='(3x,a10,"  n_trans: ",i2)',trans,d[k].n_trans
      ENDFOR 
    ENDIF
  ENDFOR
  progress,100.*float(i+1)/n
ENDFOR
progress,100.,/last

free_lun,lout

IF n_elements(problem_files) GT 1 AND NOT keyword_set(quiet) THEN BEGIN
  problem_files=problem_files[1:*]
  problem_ions=problem_ions[1:*]
  n=n_elements(problem_files)
  message,/info,/cont,'There are '+trim(n)+' files with duplicate transitions. They are listed below.'
  FOR i=0,n-1 DO BEGIN
    convertname,problem_ions[i],iz,ion
    z2element,iz-ion+1,seq,/symb
    print,'  ',problem_files[i]+'  ['+seq+' sequence]'
  ENDFOR
  message,/info,/cont,'Problem transitions are listed in '+outfile+'.'
ENDIF ELSE BEGIN 
  message,/info,/cont,'There are no files with duplicate transitions!'
ENDELSE 

END
