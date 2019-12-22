;+
; NAME:
;       READ_RRLVL
;
; PURPOSE:
;       Read the level-resolved RR rates.
;
; CATEGORY:
;       CHIANTI; data setup.
;
; CALLING SEQUENCE:
;	Result = CH_SETUP_ION( Ion_path, status )
;
; NOTES:
;       The program reads a free-format .rrlvl file, similar to the
;       .reclvl files (which are fixed format). 
;       As in the case of the .reclvl files, only the first
;       temperature array is read. 
;
; INPUTS:
;	fname:  The path name for an ion in CHIANTI format
;
; OUTPUTS:
;       A structure with the tags:
;       rate          The rate 
;       temp          The temperature (K)
;       final_level   The number identifying the final level in the
;                         recombined ion.
;       initial_level The  number identifying the initial_level in the
;                     recombining ion.
;       ref           The reference.
;
;       STATUS
;
;       A flag to tel the calling program if the .rrlvl file was
;       available (either 0 or 1)
;
; MODIFICATION HISTORY:
;       v.1, 14-Dec-2018  G. Del Zanna (GDZ)
;       v.2, 16 Jan 2019, GDZ 
;        extended to 10000 the number of lines in the files.
; 
;-

function read_rrlvl,fname, status

namerec=fname+'.rrlvl'

resultrec=findfile(expand_path(namerec))
if  resultrec(0) eq '' then begin
    status= 0
    return, -1 
 endif else status=1

    datarec=strarr(10000)
    lgt=''
    i=0
    ppp=''

    openr,1,namerec

    while strcompress(lgt,/remove_all) ne '-1' do begin 

        readf,1,lgt
;assume all temperatures are the same as in reclvl.
; read only the first one.
; free format.

        if i eq 0 then begin

         dummy= str_sep(lgt,' ',/trim)
         ind=where(trim(dummy) ne '')
         dummy= dummy[ind]
        temp=float(dummy[4:*])
        ntemp=n_elements(temp)

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

    ref=refstring(0:nref-1)
    close,1

    datarec=datarec(where(datarec ne ''))
    ntransrec=n_elements(datarec)

    final_level=intarr(ntransrec)
    initial_level=intarr(ntransrec)

    rate=dblarr(ntransrec,ntemp)
    inter=dblarr(ntemp)

aa=0 & aai=0 & aa1=0

    for i=0,ntransrec-1 do begin
        reads,datarec(i),aa,aa,aai,aa1,inter
        rate(i,*)=inter
initial_level[i]=aai
        final_level[i]=aa1
    endfor

out={rate:rate,temp:temp,$
    final_level:final_level, initial_level:initial_level,ref:ref}

return, out 


 end 
