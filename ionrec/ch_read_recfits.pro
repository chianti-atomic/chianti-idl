;        This program was developed as part of CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.


function ch_read_recfits,fits_rr,fits_dr_c,fits_dr_e


; routine for processing RR data
; read the coefficients  from the master file:

rrlines=long(file_lines(fits_rr))
rrstates=intarr(4,rrlines)
rrcoeffs=fltarr(6,rrlines)

openr,1,fits_rr

;RR RATE COEFFICIENT FITS (C)20170707 N. R. BADNELL, DEPARTMENT OF PHYSICS, UNIVERSITY OF STRATHCLYDE, GLASGOW G4 0NG, UK.;

;  Z  N  M  W      A        B        T0         T1        C        T2


; αRR(T ) = A [ (T/T0)1/2 (1+(T/T0)1/2)1-B (1+(T/T1)1/2)1+B ]-1

; B -> B + C exp(-T2/T )

str=''
bcount=0

while not eof(1) do begin

  readf,1,str
   
;  data=str_sep(trim(str),' ',/trim)
;  good=where(trim(data) ne '',nd)
;   
;  if nd eq 8 then reads,str,z,n,m,w,a,b,t0,t1 $
;  else if nd eq 10 then reads,str,z,n,m,w,a,b,t0,t1,c,t2 $
;  else begin 
;    print,'error reading the master file.. '
;    close,/all
;    return,-1
;  endelse 
   
  for istate=0,3 do rrstates[istate,bcount]=fix(strmid(str,istate*3,3))

  rrcoeffs[0,bcount]=float(strmid(str,12,11))
  rrcoeffs[1,bcount]=float(strmid(str,23,8))
  rrcoeffs[2,bcount]=float(strmid(str,31,11))
  rrcoeffs[3,bcount]=float(strmid(str,42,11))
  
  if strlen(str) gt 53 then begin
    rrcoeffs[4,bcount]=float(strmid(str,53,8))
    rrcoeffs[5,bcount]=float(strmid(str,61,11))
  endif    
   
  bcount=bcount+1
   
endwhile 

close,1

if bcount ne rrlines then message,'Error in reading RR fitting coefficients file'


; read DR data

; read the coefficients c from the master file:

drclines=long(file_lines(fits_dr_c))
drcstates=intarr(4,drclines)
drccoeffs=fltarr(9,drclines)

str=''
ccount=0

openr,2,fits_dr_c

;  Z  N  M  W      C1         C2         C3         C4         C5         C6         C7         C8         C9

while not eof(2) do begin 

  readf,2,str
   
  for istate=0,3 do drcstates[istate,ccount]=fix(strmid(str,istate*3,3))

  linelen=strlen(str)
  ncoeffs=(linelen-12)/11
  for icoeff=0,ncoeffs-1 do drccoeffs[icoeff,ccount]=float(strmid(str,icoeff*11+12,11))
   
  ccount=ccount+1
   
endwhile 

close,2

if ccount ne drclines then message,'Error in reading DR fitting coefficients C file'


drelines=long(file_lines(fits_dr_e))
drestates=intarr(4,drelines)
drecoeffs=fltarr(9,drelines)

str=''
ecount=0

openr,3,fits_dr_e

;   Z  N  M  W      E1         E2         E3         E4         E5         E6         E7         E8         E9

while not eof(3) do begin 
   
   readf,3,str
   
  for istate=0,3 do drestates[istate,ecount]=fix(strmid(str,istate*3,3))

  linelen=strlen(str)
  ncoeffs=(linelen-12)/11
  for icoeff=0,ncoeffs-1 do drecoeffs[icoeff,ecount]=float(strmid(str,icoeff*11+12,11))
   
  ecount=ecount+1
   
endwhile 

close,3

if ecount ne drelines then message,'Error in reading DR fitting coefficients E file'


if not array_equal(drcstates[0,*],drestates[0,*]) then $
  message,'Different number of C fitting coefficients than E coefficients for DR'

rec_fits={rrstates:rrstates,rrcoeffs:rrcoeffs,$
  drstates:drcstates,drcoeffs_c:drccoeffs,drcoeffs_e:drecoeffs}
  
return,rec_fits


end

