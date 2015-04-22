#!/usr/bin/perl -w 

###
### Rodger B. Voelker, U. of Oregon, Inst. of Mol. Biology
### Modified by Adam J. Struck, U. of Oregon, Dept. of Chemistry and Biochemistry
###

########### USER INPUT ###################
use strict;
use IO::File;
use PDF::API2;

### Get command line input ###
my @c_vals; # array of hashes of command line switches and values
GetARGV(\@ARGV, \@c_vals); # returns an array of hashes for each set of command line switches and values

my %switch; # hash of command line values
ProcessSwitches(\%switch, \@c_vals); # return a hash of command line switches

my $seqfile = $switch{seqfile};
my $motifs = $switch{motifs};
my $outlabel = $switch{outlabel};
###############  MAIN ################
## open file
my $seqFH = new IO::File;
$seqFH-> open("<$seqfile") or die "\nERROR: Can\'t open file \"$seqfile\"\n";

## open and read motifs file, Note pattern must be in 1st column
my $motifFH = new IO::File;
$motifFH-> open("<$motifs") or die "\nERROR: Can\'t open file \"$motifs\"\n";

my $max_motifName_len = 0;
my @motifNames;
my @motif; 

chomp($_ = <$motifFH>);  # step past header
while( <$motifFH> ) {
    my @data = split(/\t/);
    if(@data == 2) {
	push(@motif, $data[0]);

	my $motifName = $data[1];
	if(length($motifName) > $max_motifName_len) {
	    $max_motifName_len = length($motifName);
	}
	push(@motifNames, $motifName);
    }
}
$motifFH->close();

## figure out max header and sequence length
my $max_seqlen = 0;
my $max_headlen = 0;
my @headers;
my @seqs; # array of arrays of seqs
my @seqlens;
chomp($_ = <$seqFH>); # step past header
while ( <$seqFH> ) {
    my @line = split(/\t/);
    if(scalar(@line) < 2){ next; }

    my $header = $line[0];
    if(length($header) > $max_headlen){
	$max_headlen = length($header);
    }
    push(@headers, $header);
    
    my $seqlen = 0;
    my @tempseq;
    for(my $i=1; $i < scalar(@line); ++$i){ # get each seq
	my $seq = $line[$i];
	$seqlen += length($seq);
	if($seqlen > $max_seqlen) {
	    $max_seqlen = $seqlen;
	}
	push(@tempseq, $seq);
    }
    push(@seqs, \@tempseq);
    push(@seqlens, $seqlen);
}
$seqFH->close();

## open output file for pdf
my $pdfFH = new IO::File;
my $outname_pdf = $outlabel . ".pdf";
my $pdf = PDF::API2->new();
my $page = $pdf->page();
my $xpos = 0; # pos on page
my $ypos = 0;
my $font = $pdf->corefont('Helvetica', 1);	
my $fontsize = 10;
use constant in => 1/72;
my $x_left = 1/in;
my $x_right = 7.5/in; 

## establish some drawing constants
my $line_start_x = ($max_headlen*$fontsize) + 50;
my $max_line_pts = $x_right - $line_start_x;
my $pts_per_nt = $max_line_pts/$max_seqlen; # max pts left to draw the line within
my $line_center = $line_start_x + ($max_line_pts/2);
my $text = $page->text();

## my @colors = map {join "", map { sprintf "%02x", rand(255) } (0..2) } (0..63); # random color generator
my @colors = ('000000','FF0000','33FF00','0000FF', 'FFFF00');
my %colorAssignments;
assignColors(\%colorAssignments, \@motif, \@colors);

## for each string gather all starts and stops of matches, reduce to non-redundant set
# print "gene\tsegment\tmotif\tn_matches\tstarts\n";
for(my $i = scalar(@headers)-1; $i >= 0; --$i){
    if($ypos > 700){
	$page = $pdf->page(1);
	$text = $page->text();
	$ypos = 0;
    }

    ##  move to new line pos
    $ypos += 2*$fontsize;
    $xpos = $x_left;
    ## write sequence name
    $text->translate($xpos, $ypos);
    $text->font($font, $fontsize);	
    $text->text("$headers[$i]");
    $text->stroke;
    
    ## my $xstart = $line_center - ($seqlens[$i]/2);
    my $xstart = $line_start_x;
    my $xend;
    

    for(my $segment = 0; $segment < scalar(@{$seqs[$i]}); ++$segment){
	my $gfx = $page->gfx();	    

	## find matches
	my %matches; # hash of 2-element arrays, key= motif, val = non-redundant starts and stops for matches found
	my $n_matches = findMatches(\$seqs[$i][$segment], \%matches, \@motif);
	## draw line/box
	my $line_len = length($seqs[$i][$segment]) * $pts_per_nt;
	my $xend = $xstart + $line_len;
	$gfx->strokecolor('#000000');
	if($segment%2){
	    $gfx->rectxy($xstart, $ypos - ($fontsize/2), $xend, $ypos + ($fontsize/2));
	    $gfx->stroke;
	}else{
	    $gfx->move($xstart, $ypos);
	    $gfx->line($xend, $ypos);
	    $gfx->stroke;		
	}
	
	foreach my $motif (sort keys %matches) {
	    my $arrayR = $matches{$motif};
	    # print "$headers[$i]\t$segment\t$motif\t";
	    my $n_matches = $#{$arrayR} + 1;
	    # print "$n_matches\t";
	    $gfx -> fillcolor('#' . $colorAssignments{$motif}[0]);
	    $gfx -> strokecolor('#' . $colorAssignments{$motif}[0]);
	    foreach my $posR (@{$arrayR}){
		# print "@{$posR}[0],";
		my $x1 = $xstart + ($posR->[0] * $pts_per_nt);
		my $y1 = $ypos - ($fontsize/3);
		my $x2 = $xstart + ($posR->[1] * $pts_per_nt);
		my $y2 = $ypos + ($fontsize/3);
		$gfx -> rectxy($x1, $y1, $x2, $y2);
		$gfx -> fillstroke(1);
	    }
	    print "\n"
	}
	$xstart = $xend+1;
    }	
}


## Write motif key
for(my $iter = scalar(@motif)-1; $iter >= 0; --$iter){
    my $gfx = $page->gfx();
    $ypos += 2*$fontsize;
    $xpos = $x_left;

    $gfx -> fillcolor('#' . $colorAssignments{$motif[$iter]}[0]);
    $gfx -> strokecolor('#' . $colorAssignments{$motif[$iter]}[0]);
    my $x1 = $xpos;
    my $y1 = $ypos - ($fontsize/3);
    my $x2 = $xpos + 5;
    my $y2 = $ypos + ($fontsize/1.5);
    $gfx -> rectxy($x1, $y1, $x2, $y2);
    $gfx -> fillstroke(1);
    
    $text->translate($xpos + 15, $ypos);
    $text->font($font, $fontsize);	
    $text->text("$motifNames[$iter];  $motif[$iter]");
    $text->stroke;
}

$pdf->saveas("$outname_pdf" );
$pdf->end();


########  SUBROUTINES ################

sub GetARGV{
    ## does some simple error checking and assigns values to an array of hashes
    ## the resulting array then can be evaluated in a program specific fashion
	
    my $argvR = shift();
    my $switchesR = shift();
		
    if( $#{$argvR} > 100){ 
	print "\n\nFATAL ERROR: Command line overflow\n";
	PrintHelp();
	exit;
    }elsif( $#{$argvR} < 1){
	PrintHelp();
	exit;
    }elsif( ($#{$argvR}+1) % 2 ){
	print "\n\nERROR: Unmatched command line value.\n";
	PrintHelp();
	exit;
    }else{
	my $i=0;
	while( @{$argvR} ){
	    my $nextswitch = shift(@{$argvR});
	    unless( $nextswitch =~ /^-\w/ ){
		print "\nERROR: Ambiguous switch \"$nextswitch\"\n";
		PrintHelp();
		exit;
	    }
	    $switchesR->[$i]{switch} = $nextswitch;
	    $switchesR->[$i]{value} = shift(@{$argvR});
	    $i++;
	}
    }
    
    ## check for duplicate switches
    foreach my $switch (@{$switchesR}){
	my $isduplicate = -1;
	foreach my $comparison (@{$switchesR}){
	    if( $switch->{switch} eq $comparison->{switch} ){ ++$isduplicate; }
	}
	if( $isduplicate ){
	    print "\nERROR: Duplicate switch \" $switch->{switch} \" found\n";
	    PrintHelp();
	    exit;
	}
    }
}

sub PrintHelp{
	print "****************** MotifMark.pl *********************\n";
	print "Draw a diagram showing where matches to a regex are located within 1 or more sequences.\n";
	print "Note: the sequence file must have the sequence identifier located in the 1st column. All other\n";
	print "columns will be treated as sequence segments to be examined.\n";
	print "The command line switches are:     -i   name of sequence file (tab-delimited) \n";
	print "                                   -m   name of file containing motifs to match (assume header)\n";
	print "                                   -o   label for output files (\.pdf will be appended)\n";
	print "\n";

}

sub ProcessSwitches{
    my $c_valsR = pop(); #ref to command line values array
    my $switchR = pop(); #ref to switch hash
	
    ## assign default values
    $switchR->{seqfile} = '';
    $switchR->{motifs} = '';
    $switchR->{outlabel} = '';
    
    foreach my $c_hash (@{$c_valsR}){
	my $switch = $c_hash->{switch};
	my $value = $c_hash->{value};	
	
	if( $switch eq '-i' ){
	    $switchR->{seqfile} = $value;
	}elsif( $switch eq '-m' ){
	    $switchR->{motifs} = $value;
	}elsif( $switch eq '-o' ){
	    $switchR->{outlabel} = $value;
	}else{
	    print "\nERROR: Ambiguous switch \"$switch\"\n";
	    PrintHelp();
	    exit;
	}
    }
		
    ## make sure the required switches are filled and contain appropriate values
    unless( $switchR->{seqfile} ){
	print "ERROR: sequence file not specified\n";
	PrintHelp();
	exit;
    }
    unless( -e $switchR->{seqfile} ){
	print "File \" $switchR->{seqfile} \" can not be found.\n";
	PrintHelp();
	exit;
    }
    unless( $switchR->{motifs} ){
	print "ERROR: motifs file not specified\n";
	PrintHelp();
	exit;
    }
    unless( -e $switchR->{motifs} ){
	print "File \" $switchR->{seqfile} \" can not be found.\n";
	PrintHelp();
	exit;
    }
    unless( $switchR->{outlabel} ){
	print "ERROR: name of output file not specified\n";
	PrintHelp();
	exit;
    }
}

#------------- printHeaderFormatted() ------------------------
# print the header using the specified number of columns
# return the number of fields found
# yeah this was prob better done using sprintf but for now....
sub printHeaderFormatted{
    my $header = pop();
    
    my $ncolumns = 5;
    my @fields = split(/\t/, $header);
    my $nfields = scalar(@fields);

    my $nrows = int( $nfields / $ncolumns );
    if( $nfields % $ncolumns ){
	$nrows += 1;
    }
	
    my $pos;
    print "\n*************** Data Fields Available **************\n";
    my @idx;
    my @label;
    my @bL;
    my @bR;
    for( my $r = 1; $r <= $nrows; ++$r ){
	for( my $c = 0; $c < $ncolumns; ++$c ){
	    $idx[$c] = '';
	    $label[$c] = '';
	    $bL[$c] = '';
	    $bR[$c] = '';
	    $pos = $r + ($nrows * $c) -1;
	    if( $pos < $nfields ){
		$bL[$c] = '[';
		$idx[$c] = $pos+1;
		$label[$c] = $fields[$pos];
		$bR[$c] = ']';
	    }
	}
	write();
    }
#
format STDOUT = 
@<@|||@< @<<<<<<<<<<< @<@|||@< @<<<<<<<<<<< @<@|||@< @<<<<<<<<<<< @<@|||@< @<<<<<<<<<<< @<@|||@< @<<<<<<<<<<< 
$bL[0],$idx[0],$bR[0],$label[0], $bL[1],$idx[1],$bR[1],$label[1], $bL[2],$idx[2],$bR[2],$label[2], $bL[3],$idx[3],$bR[3],$label[3], $bL[4],$idx[4],$bR[4],$label[4],
.
#
	return $nfields;
}

#-------------- getInput ----------------------------
sub getInput{
    my $message = pop();
    
    print "$message";
    chomp( my $input = <STDIN> );
    return $input;
}

#-------------- getInputWithBounds ----------------------------
sub getInputWithBounds{
    my $upperbound = pop();
    my $lowerbound = pop();
    my $message = pop();
    
    my $isOK = 0;
    my $input = "";
    do{
	print "$message";
	chomp($input = <STDIN> );
	if( $input < $lowerbound or $input > $upperbound ){
	    print "**ERROR: input out of range. Acceptable range is " . $lowerbound . " to " . $upperbound . "\n";
	}else{
	    $isOK = 1;
	}
    }until( $isOK );
    return $input;
}


#------------- findMatches -------------------
sub findMatches{
    my $motifsR = pop;
    my $matchesR = pop;
    my $seqR = pop;
    my $n = 0;
    
    foreach my $m (@{$motifsR}){
	while( $$seqR =~ /($m)/g ){
	    ++$n;
	    my $mstop = pos($$seqR);
	    my $mstart = $mstop - length($1);
	    push(@{$matchesR->{$m}}, [$mstart, $mstop]);
	}
    }
    return $n;
}

#------------- assignColors -------------------
## define the colors used to mark each motif
## at some point I should rework his so that overlapping motifs are marked in a special way
## also a key would probably be nice 
sub assignColors{
    my $colorsR = pop;
    my $motif_col = pop;
    my $colorAssignmentsR = pop;
    my $n_col = 0;

    foreach my $m (@{$motif_col}) {
	my $col = @{$colorsR}[$n_col];
	push(@{$colorAssignmentsR->{$m}}, $col);
	++$n_col;
    }
}
