
FUNCTION cit_chianti_citing_papers, tot_cit=tot_cit, bcode_file=bcode_file

;+
; NAME:
;     CIT_CHIANTI_CITING_PAPERS
;
; PURPOSE:
;     Collects all the citations to the CHIANTI papers and returns a
;     structure array containing information on the papers.
;
; CATEGORY:
;     CHIANTI; citations.
;
; CALLING SEQUENCE:
;     Result = CIT_CHIANTI_CITING_PAPERS()
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     BCODE_FILE:  The name of a text file containing the
;                  bibcodes of the CHIANTI papers. The routine
;                  automatically searches for 
;                  the file from PRY's webpage, and then the
;                  CHIANTI database, so this input should only be used
;                  if you're using a different version of the file.
;	
; OUTPUTS:
;     A structure in the format returned by cit_get_ads_entry
;     containing the publication information for all of the papers
;     that cite CHIANTI.
;
; OPTIONAL OUTPUTS:
;     Tot_Cit:  The total number of citations. This is not the number
;               of elements of the output structure (which is the
;               number of *unique* citing papers). TOT_CIT is the
;               total number of citations received by the CHIANTI
;               papers, including multiple citations from the same
;               paper. 
;
; CALLS:
;     SOCK_GET, CIT_GET_ADS_ENTRY, CIT_FILL_STRINGS,
;     CIT_AFFIL_COUNTRY. 
;
; EXAMPLE:
;     IDL> s=cit_chianti_citing_papers(tot_cit=tot_cit)
;
; MODIFICATION HISTORY:
;     Ver.1, 17-Sep-2019, Peter Young
;     Ver.2, 9-Oct-2019, Peter Young
;       I had hard-coded the affil_file name - this has been removed
;       now. 
;-

;
; This reads the file that contains the list of CHIANTI paper
; bibcodes.
; By default, the copy of the routine on my website is
; read. I'll also keep a copy in the CHIANTI dbase distribution as a backup. 
;
IF n_elements(bcode_file) NE 0 THEN BEGIN
  chck=file_search(bcode_file,count=count)
  IF count EQ 0 THEN BEGIN
    print,'% CIT_CHIANTI_CITING_PAPERS: The specified BCODE_FILE does not exist. Returning...'
    return,-1
  ENDIF
  file=bcode_file
ENDIF ELSE BEGIN
  chck=have_network()
  IF chck EQ 1 THEN BEGIN
    url='http://files.pyoung.org/chianti/chianti_papers_bcodes.txt'
    sock_get,url,out_dir=getenv('IDL_TMPDIR'),local_file=local_file
    chck=file_info(local_file)
    IF chck.exists EQ 1 THEN file=local_file
  ENDIF 
  IF n_elements(file) EQ 0 THEN BEGIN 
    citdir=concat_dir(!xuvtop,'ancillary_data')
    citdir=concat_dir(citdir,'citations')
    file=concat_dir(citdir,'chianti_papers_bcodes.txt')
  ENDIF 
ENDELSE 

 
openr,lin,file,/get_lun
str1=''
WHILE eof(lin) NE 1 DO BEGIN
  readf,lin,format='(a19)',str1
  IF n_elements(bcode_list) EQ 0 THEN bcode_list=str1 ELSE bcode_list=[bcode_list,str1]
ENDWHILE 
free_lun,lin

n=n_elements(bcode_list)
biblist=''
FOR i=0,n-1 DO BEGIN
  bibs=cit_get_citing_papers(bcode_list[i])
  biblist=[biblist,bibs]
ENDFOR
biblist=biblist[1:*]

;
; Now add the extra bibcodes from chianti_extra_bibcodes.txt.
;   - these are papers that do not cite the CHIANTI papers, but
;     mention CHIANTI by name in their text and use CHIANTI for their
;     analysis. They have been manually identified by doing full-text
;     search of journals
;
citdir=concat_dir(!xuvtop,'ancillary_data')
citdir=concat_dir(citdir,'citations')
file=concat_dir(citdir,'chianti_extra_bibcodes.txt')
chck=file_info(file)
IF chck.exists EQ 1 THEN BEGIN
  print,'% CIT_CHIANTI_CITING_PAPERS: adding extra bibcodes from '+file_basename(file)+'...'
  openr,lin,file,/get_lun
  WHILE eof(lin) NE 1 DO BEGIN
    readf,lin,str1
    biblist=[biblist,trim(str1)]
  ENDWHILE 
  free_lun,lin
ENDIF 


tot_cit=n_elements(biblist)

biblist=biblist[uniq(biblist,sort(biblist))]
nb=n_elements(biblist)

print,'% CIT_CHIANTI_CITING_PAPERS: found '+trim(nb)+' unique citing papers.'

print,'% CIT_CHIANTI_CITING_PAPERS: downloading paper information from ADS...'
ads_data=cit_get_ads_entry(biblist,bad_bibcodes=bad_bibcodes)

print,'% CIT_CHIANTI_CITING_PAPERS: constructing author and article strings...'
cit_fill_strings,ads_data
;
; I've hardcoded the affiliation file location - this should be
; removed when the software has stabilized!
;
print,'% CIT_CHIANTI_CITING_PAPERS: extracting affiliation information...'
cit_affil_country,ads_data,/first

IF n_elements(bad_bibcodes) NE 0 THEN print,'% CIT_CHIANTI_CITING_PAPERS: warning - bad bibcodes found!'


return,ads_data

END
