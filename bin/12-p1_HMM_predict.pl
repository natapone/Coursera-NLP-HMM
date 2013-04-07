#!/usr/local/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use FindBin;

my $train_file_path = "$FindBin::Bin/../../h1-p/gene.train";
my $count_file_path = "$FindBin::Bin/../../h1-p/gene.counts";
#my $dev_file_path   = "$FindBin::Bin/../../h1-p/gene.dev";
my $dev_file_path   = "$FindBin::Bin/../../h1-p/gene.test";
my $assignment_name = "p1";

my $result_file_path    = "$FindBin::Bin/../../h1-p/gene_dev.$assignment_name.out";

# use result from previous script
$count_file_path = $count_file_path.".".$assignment_name;

# build library: 1-gram, 2-grams, 3-grams and word-tags

my ($lib_word_tag, $lib_n_gram, $all_n_gram) = build_lib($count_file_path);
#print Dumper($lib_word_tag), "\n";
#print Dumper($lib_n_gram), "\n";
#print Dumper($all_n_gram), "\n";

# if cant tag, get default
my $default_tag = _get_default_tag($lib_n_gram, $all_n_gram);

# Compute emission part
my $dev_strings = read_gene_dev($dev_file_path);
my $tagger_results = compute_emission($dev_strings, $lib_word_tag, $lib_n_gram, $all_n_gram, $default_tag);
#print Dumper($tagger_results) ," \n";

# write to result file
if ($assignment_name eq 'p1') {
    write_result($tagger_results, $result_file_path);
}


sub compute_emission {
    my ($dev_strings, $lib_word_tag, $lib_n_gram, $all_n_gram, $default_tag) = @_;
    
    my @tagger_results;
    
    foreach my $word (@$dev_strings) {
#        print $word, " ---- \n";
        
        my $result = {};
        $result->{'word'}   = $word;
        
        if ($word ne "") {
            $result->{'tag'}    = $default_tag;
            $result->{'score'}  = -1;
            
            # try all possible
            foreach my $tag (keys $all_n_gram->{'1'}) {
#                print "     --> ",$tag, " : ";
                
                # count(x|y) = count(y,x) / count(y)
                if ($lib_word_tag->{$tag}->{$word}->{'_count'} && $lib_n_gram->{$tag}->{'_count'}) {
                    my $scr = $lib_word_tag->{$tag}->{$word}->{'_count'} / $lib_n_gram->{$tag}->{'_count'};
                    if ($scr > $result->{'score'}) {
                        $result->{'tag'}    = $tag;
                        $result->{'score'}  = $scr;
                    }
                    
#                    print $lib_word_tag->{$tag}->{$word}->{'_count'}
#                    , " / ", 
#                    $lib_n_gram->{$tag}->{'_count'}
#                    , " = ",
#                    $lib_word_tag->{$tag}->{$word}->{'_count'} / $lib_n_gram->{$tag}->{'_count'}
#                    , "\n" ;
                } 
            }
            
        } 
        
        push(@tagger_results, $result);
    }
#    print Dumper(@tagger_results) ," \n"
    return \@tagger_results;
}

sub read_gene_dev {
    my ($dev_file_path) = @_;
    
    my @strings;
    
    # add start sequences
#    push(@strings, '*'); push(@strings, '*');
    
    my $i;
    open (MYFILE, $dev_file_path);
    while (<MYFILE>) {
        chomp;
        
        push(@strings, $_);
        
#        $i++;
#        last if ($i >= 150);
    }
    close (MYFILE);
    
    # add stop sequences
#    push(@strings, 'STOP'); push(@strings, 'STOP');
    
    return \@strings;
}

sub build_lib {
    my ($count_file_path) = @_;
    
    my $lib_n_gram      = {};
    my $lib_word_tag    = {};
    my $all_n_gram      = {}; # list of all possible result
    
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
            
            #----------
            if ($target_hash eq '$lib_n_gram' && $1 ) {
                $all_n_gram->{$1}->{$count_parts[2]} = 1;
                
            }
            #----------
            
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
    
    return ($lib_word_tag, $lib_n_gram, $all_n_gram);
}

sub write_result {
    my ($tagger_results, $result_file_path) = @_;
    
    open (OUTFILE, '>'.$result_file_path); # write to result file
    foreach my $result (@$tagger_results) {
        
#        print Dumper($result), "\n";
        
        my $out;
        if ($result->{'word'} ne "") {
            $out = $result->{'word'};
            
            if (defined($result->{'tag'})) {
                $out .= " ".$result->{'tag'};
            }
        }
        
        if (defined($out)) {
            print OUTFILE $out."\n";
        } else {
            print OUTFILE "\n";
        }
        
    }
    close (OUTFILE); 
    print "Write output to $result_file_path \n";
}


sub _get_default_tag {
    my ($lib_n_gram, $all_n_gram) = @_;
    
    my $max_count = -1;
    my $default_tag = '';
    # return most found as default
    foreach my $tag (keys $all_n_gram->{'1'}) {
        if ( $lib_n_gram->{$tag}->{'_count'} > $max_count) {
            $max_count   = $lib_n_gram->{$tag}->{'_count'};
            $default_tag = $tag;
        }
    }
    
    return $default_tag;
}

sub _simple_escape {
    my $text = shift;
    $text =~ s/'/\\'/g;
    $text =~ s/"/\\"/g;
    
    return $text;
}
