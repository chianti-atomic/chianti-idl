;+
; NAME
;
;    PROTON_DENS()
;
; EXPLANATION
;
;    Calculates the ratio of the proton density to electron density using 
;    abundance and ion balance files.
;
; INPUTS
;
;    TEMP    The logarithm (base 10) of the temperature(s) for which the 
;            ratio is required. Can be an array.
;
; OUTPUT
;
;    An array of same size as TEMP containing the proton-to-electron ratio.
;
; KEYWORDS
;
;    HYDROGEN If set then the routine computes the ratio of hydrogen to 
;             free electrons.
;
; CALLS
;
;    READ_IONEQ, READ_ABUND
;
; COMMON BLOCKS
;
;    ELEMENTS
;
; PROGRAMMING NOTES
;
;    To work out the proton/electron ratio, an ion balance and abundance 
;    file are required. These can be specified through the common block, 
;    otherwise the default files are assumed (!ioneq_file and !abund_file).
;
;    Because the ion balance data is tabulated only for logT from 4.0 to 
;    8.0, the proton/electron ratio can only be calculated for this range. 
;    Above and below these temperatures, the values at 8.0 and 4.0 are 
;    assumed, respectively.
;
;
;    I've added a check to see if the temperatures are tabulated at 0.1 
;    dex intervals (e.g., 4.0, 4.1, etc.). If they are, then a quicker 
;    algorithm is used to calculate p/e ratios. This is useful for 
;    synthetic.pro.
;
; HISTORY
;
;    Ver.1, 5-Dec-2001, Peter Young
;
;    Ver.2, 3-Dec-2001, Peter Young
;        Added /hydrogen keyword.
;
;    V. 3, 22-May-2002, GDZ: 
;                  generalize directory concatenation to work for Unix, Windows
;                  and VMS.
;
;       V.4, 06-Aug-02 GDZ
;              Changed the use of CHIANTI system variables. 
;   
; VERSION     : 4, 06-Aug-02
;
;-
FUNCTION  proton_dens, temp, hydrogen=hydrogen

common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref


IF n_elements(abund) EQ 0 THEN BEGIN 
  read_abund, !abund_file, abund,abund_ref
print,  'Using the default abundance file '+!abund_file +' to calculate the p/e ratio'
END 

IF n_elements(ioneq) EQ 0 THEN BEGIN 
  read_ioneq,!ioneq_file, ioneq_logt,ioneq,ioneq_ref
print,  'Using the default ion fraction file '+!ioneq_file+' to calculate the p/e ratio'
END 


temp=double(temp)
nt=n_elements(temp)

siz=size(ioneq)
n_ion=siz[3]

;
; tst1 is set to 1 if the temperatures are found NOT to be on (roughly) 0.1
; dex intervals.
;
tst1=0
index=0
FOR i=0,nt-1 DO BEGIN
   getmin=min( abs(temp[i]-ioneq_logt), ind)
   index=[index,ind]
   IF getmin GT 0.001 THEN tst1=1
ENDFOR
index=index[1:*]


IF tst1 EQ 0 THEN BEGIN
                                ;
   nprot=abund[0]*ioneq[index,0,1]
   nh=nprot+abund[0]*ioneq[index,0,0]
                                ;
   i_a=where(abund NE 0.,nia)
   nelec=0d0
   IF n_elements(index) EQ 1 THEN BEGIN
      FOR i=0,nia-1 DO BEGIN
         nelec=nelec+abund[i_a[i]]* $
           (transpose(reform(ioneq[index,i_a[i],*])) # findgen(n_ion))
      ENDFOR
   ENDIF ELSE BEGIN
      FOR i=0,nia-1 DO BEGIN
         nelec=nelec+abund[i_a[i]]* $
           (reform(ioneq[index,i_a[i],*]) # findgen(n_ion))
      ENDFOR
   ENDELSE

   IF keyword_set(hydrogen) THEN return,nh/nelec ELSE return,nprot/nelec
                                ;
ENDIF ELSE BEGIN
                                ;
   edens=dblarr(nt)
   pdens=dblarr(nt)
   FOR i=0,nt-1 DO BEGIN        ; loop through each temperature
      t=temp[i]
      getmin=min( abs(ioneq_logt-t), it )
      FOR j=0,siz[2]-1 DO BEGIN ; loop through each element
         ind=where(reform(ioneq[it,j,*]) NE 0.)
         IF ind[0] NE -1 THEN BEGIN
            ni=n_elements(ind)
            FOR k=0,ni-1 DO BEGIN ; loop through each ion
               y=reform(ioneq[*,j,ind[k]])
               ind2=where(y NE 0.)
                                ;
                                ; fit spline to log of ion balance data
                                ;
               y=alog10(y[ind2])
               x=ioneq_logt[ind2]
               y2=spl_init(x,y)
               yi=spl_interp(x,y,y2,t)
                                ;
                                ; for electron density, need elt. abundance, ion fraction, and 
                                ; charge of ion
                                ;
               edens[i]=edens[i]+abund[j]*10.^yi*double(ind[k])
                                ;
                                ; for proton density, just need number of protons (ionized H); if 
                                ; /hydrogen set, though, then also add number of neutral hydrogen 
                                ; atoms.
                                ;
               IF (j EQ 0) AND (ind[k] EQ 1) THEN pdens[i]=pdens[i]+abund[j]*10.^yi
               IF keyword_set(hydrogen) AND (j EQ 0) AND (ind[k] EQ 0) THEN $
                 pdens[i]=pdens[i]+abund[j]*10.^yi
            ENDFOR
         ENDIF
      ENDFOR
   ENDFOR
   return,pdens/edens
ENDELSE

END
