#!/usr/bin/perl

use warnings;
use strict;
use String::LCSS_XS qw/lcss/;
use Fatal qw/open/;
use PerlIO::gzip;

open(my $fh, "<:gzip", "/home/togoannot/plant-tm/aggregated_20120213.txt.gz");
my $lastcluster = '';
my @clstmembers;
my (%fgprnt, %nfgprnt);

while(<$fh>){
    chomp;
    my @vals = split /\t/;
    next if $vals[1] eq 'PubMed';
    if($lastcluster && ($lastcluster ne $vals[0])){
	findLCSS($lastcluster);
    }
    my @entries;
    if(index($vals[4], 'Name: Full=') > -1){
	@entries = @{ parseGI_UPT($vals[4]) };
    }elsif($vals[1] eq 'EMBL'){
	my $spos = parseEMBL($vals[4]);
	$entries[0] = ($spos > -1)?substr($vals[4], 0, $spos):$vals[4];
    }else{
	@entries = ($vals[4]);
	$entries[0] =~ s/\[[A-Z]\w+(?: \w+)*\]//;
	$entries[0] =~ s/\s*[;.]\s*$//;
    }
    @entries = grep {!/^EC=\d+(?:\.\d+)*/
			 && !/(?i-xsm:(?:un(?:(?:characteriz|nam)ed protei|know(?:n protei)?)|(?:hypothetical|expressed|putative) protei)n)/} @entries;
    for ( @entries ){
	my $fg = fingerprint($_);
	my $nfg = ngramFingerprint($_);
	$fgprnt{$fg}++;
	$nfgprnt{$nfg}++;
	print join("\t", ($vals[0], '1', $_, $fg)), "\n";
    }
    push @clstmembers, map {s/-/ /g;s/^/^/;s/  +/ /g;$_} @entries;
    $lastcluster = $vals[0];
}
close($fh);
findLCSS($lastcluster);

exit;

sub findLCSS {
    my $eid = shift;
    my %lngst;
    for(my $j = 0; $j < @clstmembers; $j++){
	for(my $k = $j+1; $k < @clstmembers; $k++){
	    my $l = lcss ($clstmembers[$j], $clstmembers[$k], 5);
	    next unless $l;
	    $l =~ s/^\s*//;
	    $l =~ s/\s*$//;
	    $lngst{$l}++;
	}
    }
    my (%fphstgrm, %fpmap);
    my (%nfphstgrm, %nfpmap);
    my @lcsm = sort {$lngst{$b}<=>$lngst{$a}} keys %lngst;
    # my @o = map {s/^\^//;$_} grep {m/^\^./ && length($_) > 5 && !m/^\^prote[a-z]*$/i && !m/^\^binding$/i} @lcsm;
    my @o = map {s/^\^//;$_} grep {m/^\^./ && length($_) > 5} @lcsm;
    if(@o){
	for ( @o ){
	    my $fg = fingerprint($_);
	    my $nfg = ngramFingerprint($_);
	    push @{ $fpmap{$fg} }, $_;
	    push @{ $nfpmap{$nfg} }, $_;
	    $fphstgrm{$fg} += $fgprnt{$fg} if $fgprnt{$fg};
	    $nfphstgrm{$nfg} += $nfgprnt{$nfg} if $nfgprnt{$nfg};
	    print join("\t", ($eid, '2f', $_, $fg)), "\n";
	    print join("\t", ($eid, '2nf', $_, $nfg)), "\n";
	}
    }else{
	print join("\t", ($eid, '2', '__N_A__')), "\n";
    }
    my $c = 0;
    for (sort {$fphstgrm{$b}<=>$fphstgrm{$a}} keys %fphstgrm){
	last if ($fphstgrm{$_} < 3 && $c > 0);
	print join("\t", ('FP', join(":", @{ $fpmap{$_} }), $fphstgrm{$_}, $fgprnt{$_}, $_)), "\n";
	$c++;
    }
    $c = 0;
    for (sort {$nfphstgrm{$b}<=>$nfphstgrm{$a}} keys %nfphstgrm){
	last if ($nfphstgrm{$_} < 3 && $c > 0);
	print join("\t", ('NFP', join(":", @{ $nfpmap{$_} }), $nfphstgrm{$_}, $nfgprnt{$_}, $_)), "\n";
	$c++;
    }
    @clstmembers = ();
    %fgprnt = ();
    %nfgprnt = ();
}

sub parseGI_UPT {
    my $w = shift;
    $w =~ s/\s*[;.]\s*$//;
    my @x = split /; /, $w;
    @x = map{ s/\w+Name: Full=//;s/Short=//;$_ } @x;
    return \@x;
}

sub parseEMBL {
    my $w = shift;
    my $curp = -1;
    my $parlv = 0;
    while($w =~ m/[();,]/g){
	my $c = substr($w, $-[0], 1);
	if($c eq ';' || $c eq ','){
	    unless($parlv){
		$curp = $-[0];
		last;
	    }
	}elsif($c eq '('){
	    $parlv++;
	}elsif($c eq ')'){
	    $parlv--;
	}
    }
    return $curp;
}

sub fingerprint {
    local $_ = shift;
    if($_){
	s/^\s*//;
	s/\s*$//;
	$_ = lc($_);
	s/[[:punct:]]|[[:cntrl:]]//g;
	return join(" ", sort split);
    }
}

sub ngramFingerprint {
    local $_ = shift;
    my $nsize = shift;
    $nsize ||= 2;
    $nsize = 2 if $nsize < 1;
    my @ngrams;
    if($_){
	$nsize--;
	$_ = lc($_);
	s/[[:punct:]]|[[:cntrl:]]|[[:space:]]//g;
	my @chars = split //;
	for(my $i = 0; $i + $nsize < @chars; $i++){
	    push @ngrams, join("", @chars[$i..($i+$nsize)]);
	}
    }
    return join(" ", sort @ngrams);
}

__END__
