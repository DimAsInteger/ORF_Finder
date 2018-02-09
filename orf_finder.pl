#!/usr/bin/perl -w
use strict;	
use warnings;	

#FOR GITHUB 02/08/18:

print "\U\t######################### ORF FINDER ######################### \n\n";
print "\UEnter fasta file:  ";
my $fastafile = <>;								
chomp $fastafile;
open(DATA, $fastafile) || die $!, print "cannot read file $fastafile \n";	
my @FASTA = <DATA>;
print @FASTA;
close DATA;

#Parse file header and save for later:

	my $header = parse_FastaSeq(@FASTA);

#Generate Forward Strand/Frames 1-3:

	my $frame1 = extract_FastaSeq(@FASTA);
	my @dna = split('', $frame1);						
	my @frames1_3  = forward_frames_1_2_3(@dna);
	
#Generate Reverse Strand/Frames 4-6:
	
	my $reverse_complement = rev_Comp($frame1);
	my @reverse_complement = split('', $reverse_complement);	
	my @frames4_6  = reverse_frames_1_2_3(@reverse_complement);

#Generate Array for Strand/Frames 1-6:
	
	my @frames1_6 = (@frames1_3, @frames4_6);
    	my @zero_placeholder = 'zero';						
    	unshift @frames1_6, @zero_placeholder; 	
    	print "\n\n\Uenter open reading frame number (1-6): ";		
    	my $frame = <>;								
    	chomp $frame;								
    	my $f = $frame;								#CM changed ($frame -1)(it was outputing seq. info. with '0' entry at STDIN to satisfy fixed array 'zero' placeholder
    	$_ = $frames1_6[$f];							

	print "\n\n\Uenter minimum ORF length (default is 50 nucleotides): ";		
    	my $entry = <>;								
    	chomp $entry;								
	my $min =50;								
	if ($entry > 50) {							
		$min = $entry							
	}									
	   	
while ($_ =~ /ATG/ig) {								

   	my $start = pos() - 3;
   		if ( /T(?:AG|GA|AA)/ig ) {					
     			my $stop = pos;
			my $lengthorf = $stop - $start;

		if ($lengthorf >= $min) { 
			my $s = substr ($_, $start, $stop - $start);
			my @seq = ();						
			@seq = ( $s =~ m/.{3}/g );				
			my $start1 = $start+1;					 
			print "\n\n$header | FRAME =  START POS = $start1, STOP POS = $stop, LEN = $lengthorf ";
			print $_." " foreach (@seq);
			print "\n\n";
			}
}	
		else {last;}
}


#SUBROUTINES:

sub parse_FastaSeq {

    	my(@File) = @_;
	my $header = '';
			
    foreach my $give_me_header (@File) {
	
	if($give_me_header =~ /^>/) {$header = $give_me_header;}			
	else {}
	}
    return $header;
}

sub extract_FastaSeq {

        my@fastaFile = @_;
        my $frame1 = '';

    foreach my $line (@fastaFile) {

        if ($line =~ /^\s*$/) {next;} #eliminate blank line
        elsif($line =~ /^\s*#/) {next;} #eliminate comment line
        elsif($line =~ /^>/) {next;}
        else {$frame1 .= $line;} # retain line and add it to the sequence string	
    }
    $frame1 =~ s/\s//g; # remove whitespace
    return ($frame1); 
}

sub rev_Comp {

	my (@rev) = @_;
	
		my $rc = join('', @rev);					
		my $reverse_complement = reverse($rc);			
		$reverse_complement =~ tr/ACGTacgt/TGCAtgca/;

	return $reverse_complement;		
}

sub forward_frames_1_2_3 {
    
    my (@framing) = @_;
            
	       $frame1 = join('', @framing); #convert array to string	 	
	    my $frame2 = substr($frame1, 1);
	    my $frame3 = substr($frame1, 2);

	    @frames1_3 = ($frame1, $frame2, $frame3);

	print "\U\nFrame 1: \n";
	print "$frame1\n\n";

	print "\U\nFrame 2: \n";
	print "$frame2\n\n";

	print "\U\nFrame 3: \n";
	print "$frame3\n\n";

        return @frames1_3;
}

sub reverse_frames_1_2_3 {
    
    my (@framing) = @_;
            
	    my $frame4 = join('', @framing); #convert array to string
            my $frame5 = substr($frame4, 1);
	    my $frame6 = substr($frame4, 2);

	    @frames4_6 = ($frame4, $frame5, $frame6);

	print "\U\nFrame 4: \n";
	print "$frame4\n\n";

	print "\U\nFrame 5: \n";
	print "$frame5\n\n";

	print "\U\nFrame 6: \n";
	print "$frame6\n\n";
            
        return @frames4_6;
}
