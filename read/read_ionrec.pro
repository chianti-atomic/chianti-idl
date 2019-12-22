
pro read_ionrec,fname,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status,$
                rec_ref,ci_ref

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
;      READ_IONREC
;
; EXPLANATION
;
;      Reads ionization and recombination total population rates from .rec and .ci files
;
; INPUTS
;
;      IZ      Atomic number of the element (i.e. Fe = 26)
;
;      ION     Ionization stage (i.e. XII = 12)
;
; OUTPUTS
;
;      REC_RATE      Total recombination rate for each level
;      ION_RATE      Collisional ionization rate for each level
;      TEMP_IONREC   Log_10 temperature (T in K)
;      LUP_REC       Upper level for the recombination rates
;      LUP_CI        Upper level for the ionization rates
;
; CALLS
;
;      ZION2FILENAME
;
; LIMITATIONS
;
;      At the moment the routine assumes that all temperatures within
;      the same file are the same. Also, that the temperatures in the
;      ionization and recombination files are the same.
;
; HISTORY
;
;      Ver.1, 23-Feb-2004, Enrico Landi (EL)
;
;      V 2, 29-Jul-2005 Giulio Del Zanna (GDZ). 
;           Modified input file names. Also rewritten the routine,
;           added some checks, and now also reads the references at
;           the end of the files.
;
;      V 3, 4-Aug-2006, Enrico Landi (EL)
;           corrected a bug in the definition of the temperature of the
;           ionization data in case the .reclvl data file is missing but
;           but the .cilvl file is available (line 222).
;
; VERSION     :   3, 4-Aug-2006
;;-


IF n_params() LT 6 THEN BEGIN
    print,'Use:  IDL> read_ionrec, fname, rec_rate, ci_rate, temp_ionrec, luprec, lupci, status,rec_ref,ci_ref'

    return
ENDIF

; Checks whether the ion/rec files exist

nameci=fname+'.cilvl'
namerec=fname+'.reclvl'

status=1.

resultci=findfile(expand_path(nameci))
resultrec=findfile(expand_path(namerec))
if resultci(0) eq '' and resultrec(0) eq '' then begin
    status=-1.
    return
endif

; reads the file *.ci

ppp=''
data=strarr(500)
lgt=''
i=0

if resultci(0) ne '' then begin

    openr,1,nameci

    while strcompress(lgt,/remove_all) ne '-1' do begin 

        readf,1,lgt

;assume all temperatures are the same !?
        if i eq 0 then begin
            ntemp=(strlen(lgt)-13)/10.           
            temp_ci=fltarr(ntemp)
            reads,lgt,aa,aa,aa,aa,temp_ci
        endif

        if  strcompress(lgt,/remove_all) ne '-1' then begin 
            readf,1,ppp
            data(i)=ppp
            i=i+1
        endif 

    endwhile  

;  get references
    refstring=strarr(500)
    nref=0
    lgt=''
    while strcompress(lgt,/remove_all) ne '-1' do begin 
        readf,1,lgt
        if strcompress(lgt,/remove_all) ne '-1' then begin 
            refstring(nref)=lgt
            nref=nref+1
        endif 
    endwhile  

    ci_ref=refstring(0:nref-1)

    close,1

    data=data(where(data ne ''))
    ntransci=n_elements(data)

    lupci=intarr(ntransci)
    ci_rate=dblarr(ntransci,ntemp)
    inter=dblarr(ntemp)
    aa1=0.

    for i=0,ntransci-1 do begin
        reads,data(i),aa,aa,aa,aa1,inter
        ci_rate(i,*)=inter
        lupci(i)=aa1
    endfor

    temp_ionrec=temp_ci

endif else begin

    ci_rate=0.0
    lupci=0.
    temp_ionrec=0.

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; reads the file *.rec    

if resultrec(0) ne '' then begin

    datarec=strarr(500)
    lgt=''
    i=0
    ppp=''

    openr,1,namerec

    while strcompress(lgt,/remove_all) ne '-1' do begin 

        readf,1,lgt
;assume all temperatures are the same !?
        if i eq 0 then begin
            ntemp=(strlen(lgt)-13)/10.           
            temp_rec=fltarr(ntemp)
            reads,lgt,aa,aa,aa,aa,temp_rec
        endif

        if  strcompress(lgt,/remove_all) ne '-1' then begin 
            readf,1,ppp
            datarec(i)=ppp
            i=i+1
        endif 

    endwhile  

;  get references
    refstring=strarr(500)
    nref=0
    lgt=''
    while strcompress(lgt,/remove_all) ne '-1' do begin 
        readf,1,lgt
        if strcompress(lgt,/remove_all) ne '-1' then begin 
            refstring(nref)=lgt
            nref=nref+1
        endif 
    endwhile  

    rec_ref=refstring(0:nref-1)
    close,1


    datarec=datarec(where(datarec ne ''))
    ntransrec=n_elements(datarec)

    luprec=intarr(ntransrec)
    rec_rate=dblarr(ntransrec,ntemp)
    inter=dblarr(ntemp)

    for i=0,ntransrec-1 do begin
        reads,datarec(i),aa,aa,aa,aa1,inter
        rec_rate(i,*)=inter
        luprec(i)=aa1
    endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 temp_ionrec=temp_rec

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endif else begin

    rec_rate=0.0
    luprec=0.
    if resultci(0) ne '' then temp_ionrec=temp_ci else temp_ionrec=0.

endelse

if n_elements(temp_ci) gt 0 and n_elements(temp_rec) gt 0 then begin 
    if max(abs(temp_ci-temp_rec)) ne 0 then begin 
        print,'% READ_IONREC: ERROR, temperatures in the two files are not the same'
        rec_rate=0.0
        luprec=0.
        temp_ionrec=0.
        status=-1
    end
end 

END
