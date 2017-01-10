#!/usr/bin/perl


#       random_dna_strings.pl
#       
#       Copyright 2011 Benjamin Tovar <scenesfromamemory4@gmail.com>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#
################################################################################
#
# NAME OF THE PROGRAM: random_dna_strings.pl
# DATE: 05/Jun/2011
# AUTHOR: Benjamin Tovar
# COMMENTS: This script takes arguments given by the user such the nucleotide frequencies of each 
# nucleotide (in a scale from 0.0 to 1.0), generates a "n" number of sequences of "n" length 
# with a FASTA header also given by the user and finally prints an output file in FASTA format
#
################################################################################

use warnings;
use strict;

print "
########################################################################

random_dna_strings.pl <- PERL SCRIPT WRITTEN BY BENJAMIN TOVAR

COMMENTS: This script takes arguments given by the user such the nucleotide frequencies of each 

nucleotide (in a scale from 0.0 to 1.0), generates a \"n\" number of sequences of \"n\" length 

with a FASTA header also given by the user and finally prints an output file in FASTA format

########################################################################\n\n";

##### seed the random stuff #####

srand(time|$$);

###### Set the number of iterations (number of random sequences to generate) ######

    print "1) Please type the number of iterations (How many random sequences do you want): 

    EXAMPLE: \"10\" \n\n";

    my $iterations = <STDIN>;
    chomp $iterations;

###### Set the length of the random DNA strings (how many nucleotides length) ######

    print "\n2) Please type the length of the random DNA strings (how many nucleotides length):

    EXAMPLE: \"50\" \n\n";

    my $length = <STDIN>;
    chomp $length;

######################### SET THE VALUE OF THE USER'S ARGUMENTS ################

# How much A% content per string:

    print "\n3) Please type the probability distribution of A content:

    REMEMBER THAT THE SUM OF THE FOUR PROBABILITIES MUST BE EQUAL TO \"1.00\"

    EXAMPLE: \"0.25\" \n\n";

    my $A_content = <STDIN>;

    print "
    ########################################################################
    # From a value of \"1.00\" as total probability, there are: ", (1-($A_content))," available
    ########################################################################\n\n";

# How much T% content per string:

    print "\n4) Please type the probability distribution of T content:

    REMEMBER THAT THE SUM OF THE FOUR PROBABILITIES MUST BE EQUAL TO \"1.00\"
    
    EXAMPLE: \"0.25\" \n\n";

    my $T_content = <STDIN>;

    print "
    ########################################################################
    # From a value of \"1.00\" as total probability, there are: ", (1-($A_content+$T_content))," available
    ########################################################################\n\n";

# How much G% content per string:

    print "\n5) Please type the probability distribution of G content:

    REMEMBER THAT THE SUM OF THE FOUR PROBABILITIES MUST BE EQUAL TO \"1.00\"

    EXAMPLE: \"0.25\" \n\n";

    my $G_content = <STDIN>;

    print "
    ########################################################################
    # From a value of \"1.00\" as total probability, there are: ", my $C_content = (1-($A_content+$T_content+$G_content))," available
    ########################################################################\n\n";

# How much C% content per string:

    print "\n6) Setting the probability distribution of C content\n\n";
    
    print $C_content,"\n";

#### Ask the user for the name of the fasta header

    print "\n7) Please, type the name of the fasta header for each sequence (is not necessary to put the >):

    EXAMPLE: \"random_seq\"\n\n";

    my $fasta_header_name =<STDIN>;

#### Ask the user for the name of output file

    print "\n8) Please, type the name of the output file:

    EXAMPLE: \"random_sequences_set.fa\"\n\n";

    my $output_file_name =<STDIN>;

####### ERASE WHITE SPACES OF THE <STDIN> INPUTS #######

    chomp ($A_content,$T_content,$G_content,$C_content,$fasta_header_name,$output_file_name);

##### Pass the values of the scalar variables to an array variable #####

    my @distribution = ($A_content,$T_content,$G_content,$C_content);

################################ RESULTS SUMMARY ###############################

	print "
------------------------------ RESULTS SUMMARY ------------------------------

SUCCESS: Here is the $iterations iterations of $length nucleotides length of 
DNA strings in FASTA format with probabilities of:

A = $A_content
T = $T_content 
C = $G_content 
G = $C_content 

EXPORTED TO FHE FILE: \"$output_file_name\" 

-----------------------------------------------------------------------------\n\n";

############################## OUTPUT FILE SETTINGS ############################

# Name of the output file

		my $output_file = "$output_file_name";

# Set the file handle "OUTPUT".

		open (OUTPUT_SEQ, ">$output_file"); 

############################## PROGRAM'S MAIN ENGINE ###########################


for(my $k=0;$k<$iterations;$k++){
	
	print OUTPUT_SEQ "\n>",$fasta_header_name,"_",($k+1),"\n";
	
		for(my $i=0;$i<$length;$i++){
	
			print OUTPUT_SEQ distribution(@distribution);
		}
}

exit;

################################################################################
############################### SUBROUTINES ####################################
################################################################################

################################################################################
# distribution
# A subroutine to generate random strings depending on the probability distribution
# of each nucleotide taken from James Tisdall's Beginning Perl for Bioinformatics
################################################################################

sub distribution{
	
	my @probability = @_;
	
	unless ($probability[0] + $probability[1] + $probability[2] + $probability[3] == 1){
		
		print "Sum of probabilites must be equal to \"1.0\"!\n";
		exit;
				}
	my $randnum = rand(1);

	if($randnum < $probability[0]) {
		return 'A';
	}elsif($randnum < $probability[0] + $probability[1]) {
		return 'T';
	}elsif($randnum < $probability[0] + $probability[1] + $probability[2]) {
		return 'C';
	}else{
		return 'G';
	}

}
