PRY's CHIANTI TOOLS
===================

This document summarizes the IDL software tools I have for formatting the CHIANTI files.


Spline fitting
==============

new_burly_ups
-------------
An updated version of the original burly_ups routine for fitting upsilons. A tutorial is available: see the file new_burly_ups_tutorial.pdf in this directory.


Converting standard file formats
================================

ADAS: chianti_adf04_v2
----------------------
Converts an ADAS adf04 format file to the CHIANTI .elvlc, .wgfa and .splups files.

IDL> chianti_adf04_v2, 'filename.adf04'


Froese Fischer
--------------
C. Froese Fischer has all of her atomic data in a standard format at her webpage:

http://www.vuse.vanderbilt.edu/~cff/mchf_collection.html

I have three routines for getting the files into CHIANTI format. Firstly you will need to copy the Froese Fischer datafile to your computer. Call this file 'ff_data.txt'.

Now read this file into an IDL structure:

IDL> process_ff_data, 'ff_data.txt', trans

Froese Fischer does not assign level indices to her data, so it is necessary to create a file that maps a CHIANTI index to the Froese Fischer filenames. An example of the level map file (which we call 'ff_level_map.txt') for Cr VIII is available. The configuration notation and level notation must be that used in the Froese Fischer file. This is the most time-consuming part of the process.

Once the level map file is created, read into IDL:

IDL> read_ff_map, 'ff_level_map.txt', levstr

Now you can create the .wgfa file with:

IDL> write_all_ff_data, trans, levstr, 'ion.wgfa'


NIST energy levels
------------------
Go to the NIST page and produce a list of energy levels in text format (not html). Cut-and-paste the file into a text file on your computer, calling it 'nist_energies.txt'. (Start the cut-and-paste from the 1st energy; don't include the titles of the columns.)

To convert this file to .elvlc format do:

IDL> nist_energies, 'nist_energies.txt', 'ion.elvlc'

By default the routine will not include columns for the theoretical energies. To add empty columns for these (i.e., all zeros), use the keyword /add_empty_theo.


Processing .wgfa files
======================

wgfa_tidy.pro
-------------
If you already have a .elvlc file correctly formatted, then you can update the wavelengths in the .wgfa file using the .elvlc energies using the routine wgfa_tidy. E.g.,

IDL> wgfa_tidy,'fe_13.wgfa','fe_13.wgfa_new',enfile='fe_13.elvlc'

You can also add the transition information to the end of each line using the keyword /ADD_INFO.


wgfa_compare.pro
----------------
This routine is useful if you have two .wgfa files for an ion and want to see how they compare. It assumes that the level indexing is the same for both files.

IDL> wgfa_compare, 'fe_13.wgfa_old', 'fe_13.wgfa_new'

The routine prints out all transitions for which there is a greater than 10% difference between the A-values. To check for transitions with a 50% difference or more, set LIMIT=0.50.

The routine also checks for any transitions that are in one file but not in another. To print out these transitions use the keyword /LIST_FILE1 for the missing transitions of FILE1 and /LIST_FILE2 for the other file.


wgfa_plot_comp.pro
------------------
This allows the A-values from two .wgfa files to be compared graphically. The plot format is that of Young (2004, A&A, 417, 785).

IDL> wgfa_plot_comp, 'fe_13.wgfa_old', 'fe_13.wgfa_new'

