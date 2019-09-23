

PRO cit_chianti_html, ads_data, outdir=outdir

;+
; NAME:
;     CIT_CHIANTI_HTML
;
; PURPOSE:
;     Creates html files containing the CHIANTI citation information
;     that are displayed on the CHIANTI web site.
;
; CATEGORY:
;     CHIANTI; citations.
;
; CALLING SEQUENCE:
;     CIT_CHIANTI_HTML
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     Outdir:  Specifies where the three output files will be written.
;	
; OUTPUTS:
;     Creates the webpages 'chianti_ADS.html' and
;     'citation_list.html', and the image 'chianti_cit_year.png'. By
;     default it places them in ~/chianti_repository/web/trunk, but
;     this can be changed by setting OUTDIR.
;
; OPTIONAL OUTPUTS:
;     Ads_Data:  The structure containing information about all the
;                papers that cite the CHIANTI papers. The format is
;                that produced by CIT_GET_ADS_ENTRY.
;
; CALLS:
;     CIT_CHIANTI_CITING_PAPERS, CIT_CHIANTI_INFO, CH_PLOT_CIT_YEAR,
;     CIT_INSTR_ALL_PAPERS 
; 
; EXAMPLE:
;     IDL> cit_chianti_html
;
; MODIFICATION HISTORY:
;     Ver.1, 17-Sep-2019, Peter Young
;-


IF n_elements(outdir) EQ 0 THEN BEGIN 
  outdir='~/chianti_repository/web/trunk'
ENDIF ELSE BEGIN
  chck=file_info(outdir)
  IF chck.exists EQ 0 THEN file_mkdir,outdir
ENDELSE 

ads_data=cit_chianti_citing_papers(tot_cit=tot_cit)

d=cit_chianti_info(ads_data)


htmlfile=concat_dir(outdir,'chianti_ADS.html')
openw,lun,htmlfile,/get_lun

n=n_elements(ads_data)
strn=strtrim(string(format='(i6)',tot_cit),2)

IF keyword_set(nrl) THEN add_str=' target="_blank"' ELSE add_str=''

printf,lun,'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN" "http://www.w3.org/TR/REC-html40/loose.dtd"> '
;
; The line below makes sure that the characters in the strings are
; printed correctly.
;
printf,lun,'<meta charset="utf-8"/>'
;
printf,lun,'<html>'
printf,lun,'<head>'
printf,lun,'<title>Citations to CHIANTI papers</title>'
printf,lun,'</head>'
printf,lun,'<body  bgcolor="#FFFFFF" vlink="#CC33CC">'
printf,lun,'<center>'
printf,lun,'<table border=0 cellpadding=0 cellspacing=0 width=800>'
printf,lun,'<tbody>'
printf,lun,'<tr><td height=30></td></tr>'
printf,lun,'<tr><td align=left>'


printf,lun,'<p><hr></p>'
printf,lun,'<h1>Citations to the CHIANTI papers</h1>'
printf,lun,'<p><hr></p>'
printf,lun,'<p>'
printf,lun,'One measure of the success of CHIANTI is through the number of citations to the CHIANTI papers in the scientific literature. The information below has been extracted from the <a href="https://ui.adsabs.harvard.edu/">ADS abstracts service</a>.'

printf,lun,'<p>Total number of citations to CHIANTI papers: <b>'+strn+'</b>*. '
printf,lun,'<p><a href=citation_list.html>Full list of papers that cite the CHIANTI papers</a>.'

w=ch_plot_cit_year(d.year,d.npapers,/buffer,/no_scale)
pngfile=concat_dir(outdir,'chianti_cit_year.png')
w.save,pngfile,resolution=96
w.close
printf,lun,'<p>Number of citations per year to the CHIANTI papers:</p>'
printf,lun,'<p align="center"><img width="600px" src="chianti_cit_year.png"></p>'



printf,lun,'<p>Summary of the major refereed journals that the citing papers have appeared in:'
printf,lun,'<p align=center><table border=1 cellpadding=3 cellspacing=0>'
printf,lun,'<tr><td><b>Journal</b></td>'+ $
     '</td><td align="center"><b>No. of papers</b>*</td></tr>'
nj=n_elements(d.abblist)
other=0
FOR i=0,nj-1 DO BEGIN
  IF d.abb_num[i] GT 0 THEN BEGIN
    printf,lun,'<tr><td>'+d.jnames[i]+'</td><td align="center">'+string(d.abb_num[i])+'</td></tr>'
  ENDIF 
ENDFOR
printf,lun,'</table></p>'

;
; Now add table with affiliations for first authors
;
n=n_elements(ads_data)
country=strarr(n)
FOR i=0,n-1 DO BEGIN
  IF ads_data[i].country.count() GT 0 THEN country[i]=ads_data[i].country[0]
ENDFOR 
k=where(country NE '',nk)
IF nk NE 0 THEN BEGIN 
  country=country[k]
  c=country[uniq(country,sort(country))]
  nc=n_elements(c)
  c_num=intarr(nc)
 ;
  printf,lun,'<p>This table shows the location of the first authors of the papers. The information is extracted from the affiliations of the authors stored in ADS. Of the '+trim(n)+' papers, '+trim(nk)+' affiliations could be extracted. There are '+trim(nc)+' countries in the list.</p>'
  printf,lun,'<p align=center><table border=1 cellpadding=3 cellspacing=0>'
  printf,lun,'<tr><td><b>Country</td><td><b>First author papers</td></tr>'
  FOR i=0,nc-1 DO BEGIN
    k=where(country EQ c[i],nk)
    c_num[i]=nk
  ENDFOR
  i=reverse(sort(c_num))
  c=c[i]
  c_num=c_num[i]
  FOR i=0,nc-1 DO BEGIN 
    printf,lun,'<tr><td>'+c[i]+'</td><td align=center>'+trim(c_num[i])+'</td></tr>'
  ENDFOR 
  printf,lun,'</table></p>'
ENDIF 


printf,lun,'<p>*- Some papers cite more than one of the CHIANTI papers, thus leading to multiple citations from a single paper. The columns above give the number of citing papers. The total number of such papers is <b>'+trim(n)+'</b>.'
printf,lun,'<p><hr>'
printf,lun,'<em>'+strtrim(systime(),2)+'</em>'


printf,lun,'</p></td></tr></tbody></table></center></body></html>'


free_lun,lun



;----------------
; WRITE OUT FILE LISTING ALL CITING PAPERS
;   - I'm using code I wrote for writing out publication lists
;     for missions and instruments.
;
data_str={instr_name: 'CHIANTI'}
htmlfile=concat_dir(outdir,'citation_list.html')
cit_instr_all_papers,ads_data,htmlfile=htmlfile,data_str=data_str, /add_year

END
