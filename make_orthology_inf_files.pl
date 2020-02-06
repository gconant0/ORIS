#!/usr/bin/perl

use strict;
use Bio::Tools::GFF;


my(@gff_files, $i, $gh_pair_file, $gh_file1, $gh_file2, $out_stub, $line, @data, $ID, %parents, @gene_chromos, $feature, $tag, $stream, @gene_starts, $gene_key,%transcript_genes,$chromo_counts, $name, $full_name, %found_trans, %found_genes, %pairs, $ka_pair_cut, $ka_dupl_cut, %dupls, $ka_val, %use_genes, $gene2_key, @genomes, %gene_transcripts, $parent, $genecnt, @aliases, @chromos, %pair_gene_hits, %dupl_gene_hits ,$tandem_file, $genome_file, %my_genes, $pairfile, $j, $last_chrom_id, $chrom_num, $ks_val, $len, $sat);

$gff_files[0]=$ARGV[0];
$gff_files[1]=$ARGV[1];
$gh_pair_file=$ARGV[2];
$gh_file1=$ARGV[3];
$gh_file2=$ARGV[4];
$ka_pair_cut=$ARGV[5];
$ka_dupl_cut=$ARGV[6];
$out_stub=$ARGV[7];

$genomes[0]="GENOME1";
$genomes[1]="GENOME2";

open(READDATA, "<" . $gh_pair_file) or die;

$line=<READDATA>;
while (!($line =~ /Gene1\s{1,6}Gene2\s{1,6}Ks/)) {
    $line=<READDATA>;
}

while ($line =<READDATA>) {
    $line=~ s/\n//;
    @data=split(/\t/, $line);
    
    if($data[2] =~ /\*/) {$sat="YES"; $data[2] =~ s/\*//;}
    else {$sat="NO";}
    
    if ($data[2] eq "NA") {$ks_val =0.0; $ka_val =0;}
    else {$ks_val=$data[2]; $ka_val =$data[3];}
    
    
    if ($data[3] <= $ka_pair_cut) {
        $found_trans{$data[0]}=1;
        $found_trans{$data[1]}=1;
        $pairs{$data[0]}{$data[1]}=$ks_val . "\t" . $ka_val . "\t" . $data[5] . "\t" . $sat;
        $pairs{$data[1]}{$data[0]}=$ks_val . "\t" . $ka_val . "\t" . $data[5] . "\t" . $sat;
    }
}
close(READDATA);
open(READDATA, "<" . $gh_file1) or die;

$line=<READDATA>;
while (!($line =~ /Gene1\s{1,6}Gene2\s{1,6}Ks/)) {
    $line=<READDATA>;
}

while ($line =<READDATA>) {
    $line=~ s/\n//;
    @data=split(/\t/, $line);
    
    if($data[2] =~ /\*/) {$sat="YES"; $data[2] =~ s/\*//;}
    else {$sat="NO";}
    
    if ($data[2] eq "IDENTICAL") {$ks_val =0.0; $ka_val=0.0;}
    else {
        if ($data[2] eq "NA"){$ks_val =0.0; $ka_val =0;}
        else {$ks_val=$data[2]; $ka_val =$data[3];}
    }
    if (exists $found_trans{$data[0]}) {
        if ($ka_val lt $ka_dupl_cut) {
            $dupls{"GENOME1"}{$data[0]}{$data[1]}=$ks_val . "\t" . $ka_val. "\t" . $data[5] . "\t" . $sat;
            $found_trans{$data[1]}=1;
        }
    }
    
    if (exists $found_trans{$data[1]}) {
        if ($ka_val lt $ka_dupl_cut) {
            $dupls{"GENOME1"}{$data[1]}{$data[0]}=$ks_val . "\t" . $ka_val. "\t" . $data[5] . "\t" . $sat;
            $found_trans{$data[0]}=1;
        }
    }
}

open(READDATA, "<" . $gh_file2) or die;

$line=<READDATA>;
while (!($line =~ /Gene1\s{1,6}Gene2\s{1,6}Ks/)) {
    $line=<READDATA>;
}

while ($line =<READDATA>) {
    $line=~ s/\n//;
    @data=split(/\t/, $line);
    
    if($data[2] =~ /\*/) {$sat="YES"; $data[2] =~ s/\*//;}
    else {$sat="NO";}
    
    if ($data[2] eq "IDENTICAL") {$ks_val =0.0; $ka_val=0.0;}
    else {
        if ($data[2] eq "NA"){$ks_val =0.0; $ka_val =0;}
        else {$ks_val=$data[2]; $ka_val =$data[3];}
    }
    
    if (exists $found_trans{$data[0]}) {
        if ($ka_val lt $ka_dupl_cut) {
            $dupls{"GENOME2"}{$data[0]}{$data[1]}=$ks_val . "\t" . $ka_val. "\t" . $data[5] . "\t" . $sat;
            $found_trans{$data[1]}=1;
        }
    }
    
    if (exists $found_trans{$data[1]}) {
        if ($ka_val lt $ka_dupl_cut) {
            $dupls{"GENOME2"}{$data[1]}{$data[0]}=$ks_val . "\t" . $ka_val. "\t" . $data[5] . "\t" . $sat;
            $found_trans{$data[0]}=1;
        }
    }
}




for($i=0; $i<2; $i++) {
    $genecnt=0;
    print "Opening GFF file ", $gff_files[$i], "\n";
    $stream =Bio::Tools::GFF->new(-file => $gff_files[$i], -gff_version => 3);
    while($feature = $stream->next_feature()) {
        # print $feature->primary_tag(), "\n"; # do something with feature
        
        if ($feature->primary_tag() eq "mRNA") {
            $name="";
            $parent="";
            foreach $tag ($feature->get_all_tags()) {
                if ($tag eq "Parent") {
                    @data = $feature->get_tag_values($tag);
                    $parent=$data[0];
                    $parent =~ s/^gene://;
                    $parent =~ s/^Gene://;
                }
                if ($tag eq "Name") {
                    @data = $feature->get_tag_values($tag);
                    $name=$data[0];
                }
                if ($tag eq "ID") {
                    @data = $feature->get_tag_values($tag);
                    $ID=$data[0];
                    $ID =~ s/^mRNA://;
                    $ID =~ s/^Transcript://;
                    $ID =~ s/^transcript://;
                }
            }
            if ($parent ne "") {
                if ($ID ne "") { $transcript_genes{$ID}=$parent;                }
            }
        }
        if ($feature->primary_tag() eq "CDS") {
            $name="";
            $parent="";
            foreach $tag ($feature->get_all_tags()) {
                if ($tag eq "Parent") {
                    @data = $feature->get_tag_values($tag);
                    $parent=$data[0];
                    $parent =~ s/^gene://;
                    $parent =~ s/^Gene://;
                    $parent =~ s/^transcript://i;
                    $parent =~ s/^Transcript://;
                }
                if ($tag eq "Name") {
                    @data = $feature->get_tag_values($tag);
                    $name=$data[0];
                }
                if ($tag eq "ID") {
                    @data = $feature->get_tag_values($tag);
                    $ID=$data[0];
                    $ID =~ s/^cds://i;
                     $ID =~ s/^CDS://i;
                }
            }
            
            if ($name eq "") {$name=$ID;}
            if ($parent ne "") {
                if ($ID ne "") {
                    $transcript_genes{$name}=$parent;
                    #print "Assiging CDS $name to $parent\n";
                }
            }
        }
    }
    $stream->close();
    print "Read transcripts from $gff_files[$i]\n";
    foreach $gene_key (keys %transcript_genes) {
        if (exists $transcript_genes{$transcript_genes{$gene_key}}) {
            $parent=$transcript_genes{$transcript_genes{$gene_key}};
            $transcript_genes{$gene_key}=$parent;
        }
        
        if (exists $found_trans{$gene_key}) {
            $found_genes{$transcript_genes{$gene_key}}=1;
        }
    }

    foreach $gene_key (keys %transcript_genes) {
        #print "Found: ", $gene_key, " => ", $transcript_genes{$gene_key}, "\n";
        $gene_transcripts{$transcript_genes{$gene_key}}{$gene_key}=1;
    }
    
    
    $stream =Bio::Tools::GFF->new(-file => $gff_files[$i], -gff_version => 3);
    $genecnt=0;
    while($feature = $stream->next_feature()) {
        # print $feature->primary_tag(), "\n"; # do something with feature
        if ($feature->primary_tag() eq "gene") {
            #print "Source tag: ", $feature->source_tag(), " at location : ", $feature->start(), " to ", $feature->end(), "\n";
            $name="";
            @aliases=();
            
            foreach $tag ($feature->get_all_tags()) {
                if ($tag eq "ID") {
                    @data = $feature->get_tag_values($tag);
                    $ID=$data[0];
                }
                if ($tag eq "Name") {
                    @data = $feature->get_tag_values($tag);
                    $name=$data[0];
                }
                if ($tag eq "Alias") {
                    @aliases=$feature->get_tag_values($tag);
                }
            }
            
            
            $full_name=$name;
            
            if (exists $found_genes{$full_name}) {
                if (!(exists $chromos[$i]{$feature->seq_id()})) {
                    $chromos[$i]{$feature->seq_id()}=$chromo_counts;
                    $chromo_counts++;
                }
                
                #print "Name: ", $name, " source tag: ", $feature->source_tag(), " at location : ", $feature->start(), " to ", $feature->end(), " ID: ", $feature->seq_id(), "\n";
                $genecnt++;
                $gene_chromos[$i]{$full_name}=$chromos[$i]{$feature->seq_id()};
                
                if ($feature->start() < $feature->end()) {  $gene_starts[$i][$chromos[$i]{$feature->seq_id()}]{$full_name}=$feature->start();}
                else { $gene_starts[$i][$chromos[$i]{$feature->seq_id()}]{$full_name}=$feature->end();}
            }
        }
        
    }
    $stream->close();
    
    $tandem_file = $out_stub . "_" . $genomes[$i] . "_tandems.txt";
    
    open(WRITEDATA, ">" . $tandem_file) or die;
    %dupl_gene_hits=();
    
    foreach $gene_key (keys %{$dupls{$genomes[$i]}}) {
        foreach $gene2_key (keys %{$dupls{$genomes[$i]}{$gene_key}}) {
            #print "Looking for pair ", $gene_key ," , ", $gene2_key, " matching genes ", $transcript_genes{$gene_key}, " amd ", $transcript_genes{$gene2_key}, "\n";
            if ((exists $gene_chromos[$i]{$transcript_genes{$gene_key}} ) && (exists $gene_chromos[$i]{$transcript_genes{$gene2_key}})) {
                if ($transcript_genes{$gene_key} ne $transcript_genes{$gene2_key} ) {
                    if (!(exists $dupl_gene_hits{$transcript_genes{$gene_key}}{$transcript_genes{$gene2_key}})) {
                        $dupl_gene_hits{$transcript_genes{$gene_key}}{$transcript_genes{$gene2_key}}=$dupls{$genomes[$i]}{$gene_key}{$gene2_key};
                        $dupl_gene_hits{$transcript_genes{$gene2_key}}{$transcript_genes{$gene_key}}=$dupls{$genomes[$i]}{$gene_key}{$gene2_key};
                    }
                    else {
                        if ($dupls{$genomes[$i]}{$gene_key}{$gene2_key} < $dupl_gene_hits{$transcript_genes{$gene_key}}{$transcript_genes{$gene2_key}}){
                                $dupl_gene_hits{$transcript_genes{$gene_key}}{$transcript_genes{$gene2_key}}=$dupls{$genomes[$i]}{$gene_key}{$gene2_key};
                                $dupl_gene_hits{$transcript_genes{$gene2_key}}{$transcript_genes{$gene_key}}=$dupls{$genomes[$i]}{$gene_key}{$gene2_key};
                        }
                    }
                }
            }
        }
    }
    
    foreach $gene_key (keys %dupl_gene_hits) {
        foreach $gene2_key (keys %{$dupl_gene_hits{$gene_key}}) {
            if ($gene_key lt $gene2_key) {
                $my_genes{$genomes[$i]}{$gene_key}=1;
                $my_genes{$genomes[$i]}{$gene2_key}=1;
                print WRITEDATA $gene_key, "\t", $gene2_key, "\t", $dupl_gene_hits{$gene_key}{$gene2_key}, "\n";
            }
        }
    }
    close(WRITEDATA);
}

$pairfile =$out_stub . "_pairs.txt";
open(WRITEDATA, ">" . $pairfile) or die;


foreach $gene_key (keys %pairs) {
    foreach $gene2_key (keys %{$pairs{$gene_key}}) {
        if ((exists $gene_chromos[0]{$transcript_genes{$gene_key}} ) && (exists $gene_chromos[1]{$transcript_genes{$gene2_key}})) {
            
            if (!(exists $pair_gene_hits{$transcript_genes{$gene_key}}{$transcript_genes{$gene2_key}})) {
                $pair_gene_hits{$transcript_genes{$gene_key}}{$transcript_genes{$gene2_key}}=$pairs{$gene_key}{$gene2_key};
               
            }
            else {
                if ($pairs{$gene_key}{$gene2_key} < $pair_gene_hits{$transcript_genes{$gene_key}}{$transcript_genes{$gene2_key}}){
                    $pair_gene_hits{$transcript_genes{$gene_key}}{$transcript_genes{$gene2_key}}=$pairs{$gene_key}{$gene2_key};
                }
            }
        }
        
    }
}

foreach $gene_key (keys %pair_gene_hits) {
    foreach $gene2_key (keys %{$pair_gene_hits{$gene_key}} ) {
        $my_genes{"GENOME1"}{$gene_key}=1;
        $my_genes{"GENOME2"}{$gene2_key}=1;
        print WRITEDATA $gene_key, "\t", $gene2_key, "\t", $pair_gene_hits{$gene_key}{$gene2_key}, "\n";
    }
}
close(WRITEDATA);

for($i=0; $i<2; $i++) {
    $genome_file= $out_stub . "_" . $genomes[$i] . "_genome.txt";
    open(WRITEDATA ,">" . $genome_file) or die;
    
    $genecnt=0;
    $chrom_num=0;
    $last_chrom_id=-1;
    for ($j=0; $j<$chromo_counts; $j++) {
        foreach $gene_key (sort {$gene_starts[$i][$j]{$a} <=> $gene_starts[$i][$j]{$b}} keys %{$gene_starts[$i][$j]}) {
            if (exists $my_genes{$genomes[$i]}{$gene_key}) {
                if ($last_chrom_id != $gene_chromos[$i]{$gene_key}) {
                    $chrom_num++;
                }
                
                print WRITEDATA $gene_key, "\t", $chrom_num, "\t", $gene_starts[$i][$j]{$gene_key}, "\n";
                $last_chrom_id=$gene_chromos[$i]{$gene_key};
                
                $genecnt++;
            }
        }
    }
    print "Wrote $genecnt genes to order file $genome_file\n";
    close(WRITEDATA);
}
