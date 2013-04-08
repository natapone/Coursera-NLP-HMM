#!/usr/local/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use FindBin;

my $train_file_path = "$FindBin::Bin/../../h1-p/gene.train";
my $count_file_path = "$FindBin::Bin/../../h1-p/gene.counts";
my $dev_file_path   = "$FindBin::Bin/../../h1-p/gene.dev";
#my $dev_file_path   = "$FindBin::Bin/../../h1-p/gene.test";
#my $assignment_name = "p1";
my $assignment_name = "p2";

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
if ($assignment_name eq 'p1') {
    my $dev_strings = read_gene_dev($dev_file_path);
    my $tagger_results = compute_emission($dev_strings, $lib_word_tag, $lib_n_gram, $all_n_gram, $default_tag);
    #print Dumper($tagger_results) ," \n";
    
    # write to result file
    write_result($tagger_results, $result_file_path);
}

if ($assignment_name eq 'p2') {
    # Compute Trigram part
    $all_n_gram = cal_trigram_TriTagger($lib_n_gram, $all_n_gram);
    print Dumper($all_n_gram);
    my $default_Tritag = _get_default_Tritag($all_n_gram);
    print "default_Tritag == $default_Tritag \n";
    
    my $dev_strings = read_gene_dev_Tritag($dev_file_path);
    my $tagger_results = cal_emission_TriTagger($dev_strings, $lib_word_tag, $lib_n_gram, $all_n_gram, $default_Tritag);
    
#    cal_emission_part
    
    
}

sub cal_emission_TriTagger {
    my ($dev_strings, $lib_word_tag, $lib_n_gram, $all_n_gram, $default_Tritag) = @_;
    
    my @tagger_results;
    
    foreach my $word_str (@$dev_strings) {
        print $word_str, " ---- \n";
        my @words = split(" ", $word_str);
#        my $result = {};
#        $result->{'word'}   = $word;
        
        # try all possible
        foreach my $tag_str (keys $all_n_gram->{'3'}) {
            print "        try - $tag_str ";
            print " = ";
            print $all_n_gram->{'3'}->{$tag_str};
            print "\n";
            
            my @tags = split(" ", $tag_str);
            
            if (@tags == @words) {
                my $gram_count = @tags;
                for (my $ii = 0; $ii < $gram_count; $ii++) {
                    print "           # ",$words[$ii] , "|", $tags[$ii];
                    print " = ";
                    
                    if ($lib_word_tag->{$tags[$ii]}->{$words[$ii]}->{'_count'} && $lib_n_gram->{$tags[$ii]}->{'_count'}) {
                        my $scr = $lib_word_tag->{$tags[$ii]}->{$words[$ii]}->{'_count'} / $lib_n_gram->{$tags[$ii]}->{'_count'};
#                        if ($scr > $result->{'score'}) {
#                            $result->{'tag'}    = $tag;
#                            $result->{'score'}  = $scr;
#                        }
                        
                        print $lib_word_tag->{$tags[$ii]}->{$words[$ii]}->{'_count'}
                        , " / ", 
                        $lib_n_gram->{$tags[$ii]}->{'_count'}
                        , " = ",
                        $lib_word_tag->{$tags[$ii]}->{$words[$ii]}->{'_count'} / $lib_n_gram->{$tags[$ii]}->{'_count'}
                        ;
                    } else {
                        print "======DATA MISSING!!!===================================";
                    }
                    print "\n";
                    
                    
                }
            } else {
                print "Error tagging count: \n", Dumper(\@tags) , " != \n", Dumper(\@words);
#                exit;
            }
            print " \n";
        }
    }
    
    
    return "";
}

sub cal_trigram_TriTagger {
    my ($lib_n_gram, $all_n_gram) = @_;
    
    # try all possible
    foreach my $tag (keys $all_n_gram->{'3'}) {
#        print "     --> ",$tag, " : \n";
        
        my @tags = split(" ", $tag);
#        print Dumper(\@tags);
        my $t_prob = cal_qML(\@tags, $lib_n_gram );
        $all_n_gram->{'3'}->{$tag} = $t_prob;
    }
    
    return $all_n_gram;
}

sub cal_qML {
    my ($tags, $lib_n_gram) = @_;
    
    my $hash_path = '';
    my $target_hash = '$lib_n_gram';
    
    # HMM Trigram part
    my $y1 = $$tags[0];
    my $y2 = $$tags[1];
    my $y3 = $$tags[2];
    
    # count(* * y1) / count(* *)
    $hash_path = str_to_hash_path("* * $y1");
    my $tr1_1 = eval $target_hash.$hash_path;
#    print "tr1_1 === $tr1_1 / ";
    
    $hash_path = str_to_hash_path("* *");
    my $tr1_2 = eval $target_hash.$hash_path;
#    print "tr1_2 === $tr1_2 \n";
    
    # count(* y1 y2) / count(* y1)
    $hash_path = str_to_hash_path("* $y1 $y2");
    my $tr2_1 = eval $target_hash.$hash_path;
#    print "tr2_1 === $tr2_1 / ";
    
    $hash_path = str_to_hash_path("* $y1");
    my $tr2_2 = eval $target_hash.$hash_path;
#    print "tr2_2 === $tr2_2 \n";
    
    # count(y1 y2 y3) / count (y1 y2)
    $hash_path = str_to_hash_path("$y1 $y2 $y3");
    my $tr3_1 = eval $target_hash.$hash_path;
#    print "tr3_1 === $tr3_1 / ";
    
    $hash_path = str_to_hash_path("$y1 $y2");
    my $tr3_2 = eval $target_hash.$hash_path;
#    print "tr3_2 === $tr3_2 \n";
    
    # count(y2 y3 STOP) / count (y2 y3)
    $hash_path = str_to_hash_path("$y2 $y3 STOP");
    my $tr4_1 = eval $target_hash.$hash_path;
#    print "tr4_1 === $tr4_1 / ";
    
    $hash_path = str_to_hash_path("$y2 $y3");
    my $tr4_2 = eval $target_hash.$hash_path;
#    print "tr4_2 === $tr4_2 \n";
    
    # validate 
    my $t_prob = 0;
    if (
        defined($tr1_1) && defined($tr1_2) && 
        defined($tr2_1) && defined($tr2_2) && 
        defined($tr3_1) && defined($tr3_2) && 
        defined($tr4_1) && defined($tr4_2)
    )  {
        $t_prob =    ($tr1_1 / $tr1_2) * 
                    ($tr2_1 / $tr2_2) * 
                    ($tr3_1 / $tr3_2) * 
                    ($tr4_1 / $tr4_2)
    }
    
    # assign prob for each gram
    return $t_prob;
}

sub str_to_hash_path {
    my ($str) = @_;
    
    my $path = "";
    my @parts = split(' ', $str);
    
    foreach my $part (@parts) {
#        print $arr ," -- ";
        $path .= "->{'$part'}";
    }
    return $path."->{'_count'}";
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

sub read_gene_dev_Tritag {
    my ($dev_file_path) = @_;
    
    my @strings;
    
    # add start sequences
    push(@strings, '*'); push(@strings, '*');
    
    my $i;
    open (MYFILE, $dev_file_path);
    while (<MYFILE>) {
        chomp;
        $_ = '\n' if ($_ eq ''); # prevent error splitting
        push(@strings, $_);
        
        $i++;
        last if ($i >= 30);
    }
    close (MYFILE);
    
#    print Dumper(\@strings);
    
    # return 3-gram
    my $arr_count = @strings;
    my @trigrams;
    for (my $ii = 0; $ii < $arr_count; $ii++) {
        if ( $ii-2 >=0 ) {
            push(@trigrams, $strings[$ii-2]." ".$strings[$ii-1]." ".$strings[$ii]);
        }
    }
    
#print Dumper(\@trigrams);
    
    
    
    return \@trigrams;
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
#        last if ($i >= 100);
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
            
            #-Build all possible grams combination---------
            if ($target_hash eq '$lib_n_gram' && $1 ) {
                # remove start, stop sequences from trigram possible results
                if ($count_parts[2] !~ /\*/ && $count_parts[2] !~ /STOP/) {
                    $all_n_gram->{$1}->{$count_parts[2]} = 1;
                }
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

#for Trigram tagger
sub _get_default_Tritag {
    my ($all_n_gram) = @_;
    
    return ( (sort { $all_n_gram->{'3'}->{$b} <=> $all_n_gram->{'3'}->{$a} } keys $all_n_gram->{'3'})[0] );
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
