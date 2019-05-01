#!/usr/bin/env perl
$latex              = 'uplatex -synctex=1 -halt-on-error';
$latex_silent       = 'uplatex -synctex=1 -halt-on-error -interaction=batchmode';
$bibtex             = 'pbibtex';
$dvipdf             = 'dvipdfmx %O -o %D %S';
$makeindex          = 'mendex %O -o %D %S';
$max_repeat         = 5;
$pdf_mode           = 3; # generates pdf via dvipdfmx

# Prevent latexmk from removing PDF after typeset.
$pvc_view_file_via_temporary = 0;

# Use Evince editor
$pdf_previewer    = "evince %O %S";

