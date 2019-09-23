
This document explains how to use software for processing the Froese
Fischer radiative data into CHIANTI format.

The data are available at:

http://nlte.nist.gov/MCHF/

and the data-sets you need to look for are 'Transition rates'. The
data are given in ascii files that (mostly) have the same format.

Generally the user will have created a .elvlc for the ion using the
level indexing from the collision data file and observed energies from
NIST. The Froese Fischer data should only be processed once the .elvlc file
is finalized.

A key step (no. 4) is to manually edit the level mapping from the Froese
Fischer levels to the .elvlc levels. 


The procedure is:

1. Download the entire text file to your working directory.


2. Run the routine ff_read_data to read the data file into an IDL
structure:

IDL> ff_read_data, filename, trans


3. Extract the list of levels from this structure into a text file:

IDL> ff_make_level_list, trans, outfile=outfile

The level indices (column 3) are ordered according to level energy of
the Froese Fischer model.


4. Manually edit the level indexing for this file by comparing with
the CHIANTI elvlc file. This creates a new file 'outfile_new'. 


5. Read the new level indexing file into IDL:

IDL> ff_read_levels, outfile_new, levstr


6. Create a CHIANTI format .wgfa file:

IDL> ff_write_wgfa, trans, levstr, 'ion.wgfa'
