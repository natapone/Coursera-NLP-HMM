#!/usr/local/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use FindBin;

my $train_file_path = "$FindBin::Bin/../../h1-p/gene.train";
my $count_file_path = "$FindBin::Bin/../../h1-p/gene.counts";
my $assignment_name = "p1";

# use result from previous script
$count_file_path = $count_file_path.".".$assignment_name;

# build library: 1-gram, 2-grams, 3-grams and word-tags

my ($lib_word_tag, $lib_n_gram) = build_lib($count_file_path);
print Dumper($lib_word_tag), "\n";
print Dumper($lib_n_gram), "\n";






sub build_lib {
    my ($count_file_path) = @_;
    
    my $lib_n_gram = {};
    my $lib_word_tag = {};
    
    my $i;
    open (MYFILE, $count_file_path);
    while (<MYFILE>) {
        chomp;
        my @count_words = split(/ /, $_);
        
        my $validate_pass = 0;
        my $target_hash   = '';
        # All-grams
        if ($count_words[1] =~ /(\d)-GRAM/) { # validate result
            
            if (@count_words < ($1 + 2)) {
                print "DIE @ $1-grams => $_ \n";
                exit;
            }
            
            # already pass validation
            $validate_pass = 1;
            $target_hash   = '$lib_n_gram';
            
        } elsif ($count_words[1] =~ /WORDTAG/) {
            if (@count_words < (4)) {
                print "DIE @ WORDTAG => $_ \n";
                exit;
            }
            
            # already pass validation
            $validate_pass = 1;
            $target_hash   = '$lib_word_tag';
            
        }
        
        if ($validate_pass && $target_hash) {
            my @count_parts = split(/ /, $_, 3);
            my @word_parts  = split(/ /, $count_parts[2]);
            
            my $hash_path = '';
            foreach(@word_parts) {
                $hash_path .= "->{'"._simple_escape($_)."'}";
            }
            #            print "hash_path = $hash_path \n";
            
            my $eval_string = $target_hash.$hash_path."->{'_count'} = ".$count_parts[0].";";
#            print "***eval_string === ",$eval_string, "\n";
            eval $eval_string ; die $@ if $@;
        }
        
#        $i++;
#        last if ($i >= 1000);
    }
    close (MYFILE); 
    
#    print Dumper($lib_word_tag), "\n";
#    print Dumper($lib_n_gram), "\n";
    
    return ($lib_word_tag, $lib_n_gram);
}

sub _simple_escape {
    my $text = shift;
    $text =~ s/'/\\'/g;
    $text =~ s/"/\\"/g;
    
    return $text;
}
