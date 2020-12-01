

function cit_chianti_info, ads_data, quiet=quiet

;+
; NAME:
;    CIT_CHIANTI_INFO
;
; PURPOSE:
;    Extracts information about papers citing papers that is used on
;    the CHIANTI website.
;
; CATEGORY:
;    CHIANTI; citations.
;
; CALLING SEQUENCE:
;    Result = CIT_CHIANTI_INFO( Ads_Data )
;
; INPUTS:
;    Ads_Data:  A structure in the format returned by
;               CIT_GET_ADS_ENTRY containing information about the
;               papers that cite CHIANTI.
;
; KEYWORD PARAMETERS:
;     QUIET:  If set, then no information is printed to the screen.
;
; OUTPUTS:
;     A structure is returned with the following tags:
;      year - array with a list of years (integers)
;      npapers - same size as year containing citations for each year
;      abblist - list of bibcode journal abbreviations
;      abb_num - same size as abblist containing number of papers from
;                the journals.
;      jnames - same size as abblit with names of journals (shorthand
;               format).
;
;     The routine also prints to the screen a list of journals with
;     the citation numbers, that are not part of the journal
;     list. Either the numbers are too low, or they aren't
;     journals (e.g., preprints, conference proceedings).
;
; EXAMPLE:
;     IDL> s=cit_chianti_citing_papers()
;     IDL> d=cit_chianti_info(s)
;
; MODIFICATION HISTORY:
;     Ver.1, 17-Sep-2019, Peter Young
;-


yr=fix(ads_data.year)

t=systime(/jul,/utc)
caldat,t,dd,mm,yy

h=histogram(yr,min=1995,max=yy,loc=loc,bin=1)

;
; Extract the 5 digit journal IDs from the bibcodes.
;
abb=strmid(ads_data.bibcode,4,5)

;
; Create list of the journals to be checked. ABBLIST contains the
; 5-digit codes, and JNAMES contains the names that will be printed to
; the webpage. These arrays should match!
; 
abblist=['Natur','Sci..','NatCo','NatAs','SciA.','ApJ..','ApJS.','AJ...','A&A..','MNRAS', $
         'SoPh.','ADNDT','AdSpR','PASJ.','SSRv.','AstL.','JGRA.','JPhB.','PhyS.','Icar.','JGR..']
jnames=['Nature','Science','Nat. Comm.','Nat. Ast.','Sci. Adv.','ApJ','ApJS','AJ','A&A','MNRAS', $
        'Sol. Phys.','ADNDT','Adv. Sp. Res.','PASJ','Sp. Sci. Rev.','Ast. Lett.','JGR','J. Phys. B','Phys. Scripta','Icarus','JGR' ]

m=n_elements(ads_data)
flag=bytarr(m)

n=n_elements(abblist)
num=intarr(n)
FOR i=0,n-1 DO BEGIN
  k=where(abb EQ abblist[i],nk)
  num[i]=nk
  flag[k]=1b
ENDFOR

k=where(flag EQ 0)
a=ads_data[k]
abb=abb[k]

i=uniq(abb,sort(abb))
n=n_elements(i)
FOR j=0,n-1 DO BEGIN
  k=where(abb[i[j]] EQ abb,nk)
  IF nk GT 5 AND NOT keyword_set(quiet) THEN print,'Missing abbrev: ',abb[i[j]],nk
ENDFOR 

output={ year: loc, $
         npapers: h, $
         abblist: abblist, $
         abb_num: num, $
         jnames: jnames}

return,output



END
