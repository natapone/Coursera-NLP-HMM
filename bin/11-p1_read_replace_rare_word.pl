#!/usr/local/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use FindBin;


my $train_file_path = "$FindBin::Bin/../../h1-p/gene.train";
my $count_file_path = "$FindBin::Bin/../../h1-p/gene.counts";
my $assignment_name = "p1";


print "count file name = $count_file_path \n";

print "Reading list of RARE words ..... \n";
# get list of RARE words
my $rare_words = get_rare_words($count_file_path, 5);
#my @rare_count = keys $rare_words;
print "RARE words count = ", scalar keys $rare_words, "\n";
#print Dumper($rare_words), "\n";


print "Replacing RARE words on training set ..... \n";
# replace RARE words on training set
replace_training_set ($train_file_path, $rare_words,$assignment_name);

sub replace_training_set {
    my ($train_file_path, $rare_words, $assignment_name) = @_;
    
    my $i;
    open (MYFILE, $train_file_path);
    open (OUTFILE, '>'.$train_file_path.".".$assignment_name); # new training file
    while (<MYFILE>) {
        chomp;
        my @count_words = split(/ /, $_);
        
        if (@count_words > 2) {
            print "Sucka!!! ==> $_ \n";
            exit;
        }
        
        my $word = $count_words[0];
        my $tag  = $count_words[1];
        
        my $out  = '';
        if ($word) {
            if ($rare_words->{$word}) {
                $out = "_RARE_"." ".$tag;
            } else {
                $out = $_ if ($_);
            }
        }
        
        print OUTFILE "$out\n";
        
#        $i++;
#        last if ($i >= 1500);
    }
    
    close (MYFILE);
    close (OUTFILE); 
}

sub get_rare_words {
    my ($count_file_path, $max_f) = @_;
    
    my $rare_words = {};
    my $i;
    open (MYFILE, $count_file_path);
    while (<MYFILE>) {
        chomp;
        my @count_words = split(/ /, $_);
#        print Dumper(@count_words), "\n";
        
        my $str_out;
        # use only WORDTAG
        if ($count_words[1] eq 'WORDTAG') {
            
            # if frequency < max_f --> rare words
            if ($count_words[0] < $max_f) {
#                print $count_words[0] , " of ", $max_f, " === ", $count_words[3], "\n";
                $rare_words->{$count_words[3]} = 1;
                
            }
        }
        
#        $i++;
#        last if ($i >= 1000);
    }
    close (MYFILE); 
    
    return $rare_words;
}
