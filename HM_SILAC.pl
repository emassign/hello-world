#!/usr/bin/perl 
use strict;
use LIB::modifications;
use LIB::maxquant;


# INPUT
my $usage  = 'HM_SILAC.pl <config>';
	die $usage if (!defined $ARGV[0]);
my $RT_threshold;
my $mass_threshold;
my $Score_Min;
my $Delta_Min;
my $LcPrb_Min;
my $Mod_Window;
my $database;
my (%Non_Redundant, %Correction, %Correction_data, %Results);
my (@files);

########### READ CONFIGURATION FILE ##############
open CONFIG, $ARGV[0] or die "Can't open $ARGV[0]: $!";
while (my $line = <CONFIG>)
	{
	chomp ($line);
	my @in = split ('\t', $line);
	$RT_threshold 	=	$in[1] if ($in[0] =~ /^Retention/);
	$mass_threshold =	$in[1] if ($in[0] =~ /^Mass/);
	$Score_Min		=	$in[1] if ($in[0] =~ /^Minimum Score/);
	$Delta_Min		=	$in[1] if ($in[0] =~ /^Minimum Delta/);
	$LcPrb_Min		=	$in[1] if ($in[0] =~ /^Minimum Localization/);
	$database		=	$in[1] if ($in[0] =~ /^FASTA/);
	$Mod_Window		=	$in[1] if ($in[0] =~ /^Modification/);
	if ($in[0] =~ /^Input/)
		{
		shift (@in);
		@files = @in;
		}
	}
close CONFIG;
die 'Error in config file' if (!defined $database || scalar(@files)==0);
$RT_threshold	=	0.5 unless (defined $RT_threshold);
$mass_threshold	=	0.002 unless (defined $mass_threshold);
$Score_Min	=	25 unless (defined $Score_Min);
$Delta_Min	=	20 unless (defined $Delta_Min);
$LcPrb_Min	=	0.75 unless (defined $LcPrb_Min);
$Mod_Window	=	31 unless (defined $Mod_Window);

my $parameters = $RT_threshold."min_".$mass_threshold."Da";
my $dir		= 'HM_SILAC_input/';
my $fasta	= 'Database/'.$database.".fasta";
my %fasta_database = %{&maxquant::read_and_clean_fasta(\$fasta)};
# OUTPUT
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$mon++;
my $date		= sprintf("%02d", $year % 100).'_'.$mon.'_'.$mday;
my $header_red 	= maxquant::print_header_redundant;
my $header_nr 	= maxquant::print_header_nr;
my $output_folder = "HM_SILAC_output/";


########## DATA INPUT ############

print "\nBEGIN\n\n";
# CYCLE THROUGH PROCESSED FILE DIRECTORY
#READ INPUT FILES (msms;evidence) 
foreach my $file (@files)
{
	my $subfolder = $file."_".$date."_".$hour.$min;
	
	mkdir ($output_folder.$subfolder);
	
	print "stampo ".$file."_combined\n";
	my (%msms_info,%msms_columns,%evidence_columns,%mod_peptide_summary);
	
		# OPEN/PREPARE RUNDUNANT FILE FOR OUTPUT
		open OUT_REDUNDANT, '>'.$output_folder.$subfolder.'/peptides_summary.txt';
		print OUT_REDUNDANT "$header_red\n";

	
	my $msms = $dir.$file."/msms.txt";
	open MSMS, $msms or die "Could not open $msms --> $!";
	while (<MSMS>)
	{
		chomp($_);
		my (@columns) = split(/\t/,$_);
		
		if ($_ =~ /^Raw/)
		{
			%{$msms_columns{$file}} = %{&maxquant::columns(\@columns)};			
		}
		else
		{
				next if ($columns[$msms_columns{$file}{"Reverse"}] =~ /\+/);	#FILTERING
		
				my ($peak_list,$peak_list_corrected);
				
				my (@masses) 		= split(/\;/,$columns[$msms_columns{$file}{"Masses"}]);
				my (@errors) 		= split(/\;/,$columns[$msms_columns{$file}{"Mass Deviations [Da]"}]);
				my (@intensities) 	= split(/\;/,$columns[$msms_columns{$file}{"Intensities"}]);
		
				for (my $i=0;$i<@masses;$i++)
				{
					$peak_list .= $masses[$i]."_".$intensities[$i].";";
					$peak_list_corrected .= ($masses[$i]+$errors[$i])."\s".$intensities[$i].";";
				}
				
				$msms_info{$columns[$msms_columns{$file}{"id"}]}{'Peak List'} = $peak_list;
				$msms_info{$columns[$msms_columns{$file}{"id"}]}{'Peak List Corrected'} = $peak_list_corrected;
				$msms_info{$columns[$msms_columns{$file}{"id"}]}{'Peak Coverage'} = $columns[$msms_columns{$file}{"Peak coverage"}];
				
		}
	};
	close MSMS;

	open EVIDENCE, $dir.$file."/evidence.txt" or die "Could not open ".$file;
	while (<EVIDENCE>)
	{
		chomp($_);
		my (@columns) = split(/\t/,$_);
	
		if ($_ =~ /^Sequence/)
		{
			$_ =~ s/Gene names/Gene Names/;
			$_ =~ s/scan number/Scan Number/;
			$_ =~ s/Protein names/Protein Names/;	#this is needed to work on newer versions of MaxQuant as the header of the output tables differs
			@columns = split(/\t/,$_);
			%{$evidence_columns{$file}} = %{&maxquant::columns(\@columns)};
			die 'Check evidence.txt columns names' unless (exists $evidence_columns{$file}{"Modified sequence"} && exists $evidence_columns{$file}{"Gene Names"});
		}
		else
		{
			next if ($columns[$evidence_columns{$file}{"Potential contaminant"}] =~ /\+/);		#PEPTIDE FILTERING
			next if ($columns[$evidence_columns{$file}{"Reverse"}] =~ /\+/);
			my $experiment;	
			$experiment = $columns[$evidence_columns{$file}{"Experiment"}];
			
			#PEPTIDE FILTERING 2
			next if ($columns[$evidence_columns{$file}{"OxMet4"}] != 0 && $columns[$evidence_columns{$file}{"Methyl (KR)"}] == 0 && $columns[$evidence_columns{$file}{"Dimethyl (KR)"}] == 0 && $columns[$evidence_columns{$file}{"Trimethyl (K)"}] == 0 && $columns[$evidence_columns{$file}{"Methyl4 (KR)"}] == 0 && $columns[$evidence_columns{$file}{"Dimethyl4 (KR)"}] == 0 && $columns[$evidence_columns{$file}{"Trimethyl4 (K)"}] == 0 && $columns[$evidence_columns{$file}{"Methyl (KRDEQH)"}] == 0 && $columns[$evidence_columns{$file}{"Methyl4 (KRDEQH)"}] == 0);
			next if ($columns[$evidence_columns{$file}{"Score"}] < $Score_Min);	
			next if ($columns[$evidence_columns{$file}{"Delta score"}] < $Delta_Min);
	
			my ($Methyl_all,$Methyl4_all,$Methyl,$Dimethyl,$Trimethyl,$Methyl4,$Dimethyl4,$Trimethyl4,$Met4,$Oxidation,$Met4ox) = (0,0,0,0,0,0,0,0,0,0,0);
			my ($Methyl_all_prb,$Methyl4_all_prb,$Methyl_prb,$Dimethyl_prb,$Trimethyl_prb,$Methyl4_prb,$Dimethyl4_prb,$Trimethyl4_prb,$Met4_prb,$Oxidation_prb,$Met4ox_prb) = ('','','','','','','','','','','');
				
			$Methyl_all		= $columns[$evidence_columns{$file}{"Methyl (KRDEQH)"}] if (exists($evidence_columns{$file}{"Methyl (KRDEQH)"}) && $columns[$evidence_columns{$file}{"Methyl (KRDEQH)"}] =~ /^\d+$/);
			$Methyl4_all	= $columns[$evidence_columns{$file}{"Methyl4 (KRDEQH)"}] if (exists($evidence_columns{$file}{"Methyl4 (KRDEQH)"}) && $columns[$evidence_columns{$file}{"Methyl4 (KRDEQH)"}] =~ /^\d+$/);
			$Methyl			= $columns[$evidence_columns{$file}{"Methyl (KR)"}] if (exists($evidence_columns{$file}{"Methyl (KR)"}) && $columns[$evidence_columns{$file}{"Methyl (KR)"}] =~ /^\d+$/);
			$Dimethyl		= $columns[$evidence_columns{$file}{"Dimethyl (KR)"}] if (exists($evidence_columns{$file}{"Dimethyl (KR)"}) && $columns[$evidence_columns{$file}{"Dimethyl (KR)"}] =~ /^\d+$/);
			$Trimethyl		= $columns[$evidence_columns{$file}{"Trimethyl (K)"}] if (exists($evidence_columns{$file}{"Trimethyl (K)"}) && $columns[$evidence_columns{$file}{"Trimethyl (K)"}] =~ /^\d+$/);
			$Methyl4		= $columns[$evidence_columns{$file}{"Methyl4 (KR)"}] if (exists($evidence_columns{$file}{"Methyl4 (KR)"}) && $columns[$evidence_columns{$file}{"Methyl4 (KR)"}] =~ /^\d+$/);
			$Dimethyl4		= $columns[$evidence_columns{$file}{"Dimethyl4 (KR)"}] if (exists($evidence_columns{$file}{"Dimethyl4 (KR)"}) && $columns[$evidence_columns{$file}{"Dimethyl4 (KR)"}] =~ /^\d+$/);
			$Trimethyl4		= $columns[$evidence_columns{$file}{"Trimethyl4 (K)"}] if (exists($evidence_columns{$file}{"Trimethyl4 (K)"}) && $columns[$evidence_columns{$file}{"Trimethyl4 (K)"}] =~ /^\d+$/);
			$Met4			= $columns[$evidence_columns{$file}{"Met4"}] if (exists($evidence_columns{$file}{"Met4"}) && $columns[$evidence_columns{$file}{"Met4"}] =~ /^\d+$/);
			$Met4ox			= $columns[$evidence_columns{$file}{"OxMet4"}] if (exists($evidence_columns{$file}{"OxMet4"}) && $columns[$evidence_columns{$file}{"OxMet4"}] =~ /^\d+$/);
			$Oxidation		= $columns[$evidence_columns{$file}{"Oxidation (M)"}] if (exists($evidence_columns{$file}{"Oxidation (M)"}) && $columns[$evidence_columns{$file}{"Oxidation (M)"}] =~ /^\d+$/);
			
			$Methyl_all_prb	= $columns[$evidence_columns{$file}{"Methyl (KRDEQH) Probabilities"}] if (exists($evidence_columns{$file}{"Methyl (KRDEQH) Probabilities"}) && $columns[$evidence_columns{$file}{"Methyl (KRDEQH) Probabilities"}] !~ /^$/);
			$Methyl4_all_prb= $columns[$evidence_columns{$file}{"Methyl4 (KRDEQH) Probabilities"}] if (exists($evidence_columns{$file}{"Methyl4 (KRDEQH) Probabilities"}) && $columns[$evidence_columns{$file}{"Methyl4 (KRDEQH) Probabilities"}] !~ /^$/);
			$Methyl_prb		= $columns[$evidence_columns{$file}{"Methyl (KR) Probabilities"}] if (exists($evidence_columns{$file}{"Methyl (KR) Probabilities"}) && $columns[$evidence_columns{$file}{"Methyl (KR) Probabilities"}] !~ /^$/);
			$Dimethyl_prb	= $columns[$evidence_columns{$file}{"Dimethyl (KR) Probabilities"}] if (exists($evidence_columns{$file}{"Dimethyl (KR) Probabilities"}) && $columns[$evidence_columns{$file}{"Dimethyl (KR) Probabilities"}] !~ /^$/);
			$Trimethyl_prb	= $columns[$evidence_columns{$file}{"Trimethyl (K) Probabilities"}] if (exists($evidence_columns{$file}{"Trimethyl (K) Probabilities"}) && $columns[$evidence_columns{$file}{"Trimethyl (K) Probabilities"}] !~ /^$/);
			$Methyl4_prb	= $columns[$evidence_columns{$file}{"Methyl4 (KR) Probabilities"}] if (exists($evidence_columns{$file}{"Methyl4 (KR) Probabilities"}) && $columns[$evidence_columns{$file}{"Methyl4 (KR) Probabilities"}] !~ /^$/);
			$Dimethyl4_prb	= $columns[$evidence_columns{$file}{"Dimethyl4 (KR) Probabilities"}] if (exists($evidence_columns{$file}{"Dimethyl4 (KR) Probabilities"}) && $columns[$evidence_columns{$file}{"Dimethyl4 (KR) Probabilities"}] !~ /^$/);
			$Trimethyl4_prb	= $columns[$evidence_columns{$file}{"Trimethyl4 (K) Probabilities"}] if (exists($evidence_columns{$file}{"Trimethyl4 (K) Probabilities"}) && $columns[$evidence_columns{$file}{"Trimethyl4 (K) Probabilities"}] !~ /^$/);
			$Met4_prb		= $columns[$evidence_columns{$file}{"Met4 Probabilities"}] if (exists($evidence_columns{$file}{"Met4 Probabilities"}) && $columns[$evidence_columns{$file}{"Met4 Probabilities"}] !~ /^$/);
			$Met4ox_prb		= $columns[$evidence_columns{$file}{"OxMet4 Probabilities"}] if (exists($evidence_columns{$file}{"OxMet4 Probabilities"}) && $columns[$evidence_columns{$file}{"OxMet4 Probabilities"}] !~ /^$/);
			$Oxidation_prb	= $columns[$evidence_columns{$file}{"Oxidation (M) Probabilities"}] if (exists($evidence_columns{$file}{"Oxidation (M) Probabilities"}) && $columns[$evidence_columns{$file}{"Oxidation (M) Probabilities"}] !~ /^$/);
			
			next if (($Methyl_all > 0 || $Methyl > 0 || $Dimethyl > 0 || $Trimethyl > 0) && ($Methyl4_all > 0 || $Methyl4 > 0 || $Dimethyl4 > 0 || $Trimethyl4 > 0 || $Met4 > 0 || $Met4ox > 0));
			next if ($columns[$evidence_columns{$file}{"Sequence"}] =~ /M/g && ($Methyl4_all > 0 || $Methyl4 > 0 || $Dimethyl4 > 0 || $Trimethyl4 > 0) && ($Met4 == 0 && $Met4ox == 0));
			next if ($columns[$evidence_columns{$file}{"Sequence"}] =~ /M/g && ($Methyl_all > 0 || $Methyl > 0 || $Dimethyl > 0 || $Trimethyl > 0) && ($Met4 > 0 || $Met4ox == 0));		
			next if ($columns[$evidence_columns{$file}{"Charge"}] < 2 || !defined($columns[$evidence_columns{$file}{"Charge"}]));
										
			my $Modified_Sequence = $columns[$evidence_columns{$file}{"Modified sequence"}] if exists $evidence_columns{$file}{"Modified sequence"};
			my @proteins = split(";", $columns[$evidence_columns{$file}{"Proteins"}]);
			my $leading_proteins = @proteins[0];
	
			my ($localised) = &modifications::process_localisation(
																	\$Modified_Sequence,
																	\$leading_proteins,
																	\$Methyl_all_prb,
																	\$Methyl4_all_prb,
																	\$Methyl_prb,
																	\$Dimethyl_prb,
																	\$Trimethyl_prb,
																	\$Methyl4_prb,
																	\$Dimethyl4_prb,
																	\$Trimethyl4_prb,
																	\$Oxidation_prb,
																	\$Met4_prb,
																	\$Met4ox_prb,
																	\$LcPrb_Min,
																	\%fasta_database,
																	\$Mod_Window
																);							#serve contare le modfiche su ogni peptide basandosi sulla probabilità con cui sono assegnate le modfiche
			
#PEPTIDE FILTERING 3		
			my $Kme_count = keys %{$$localised{'GOOD'}{'K'}{'me'}};
			my $Rme_count = keys %{$$localised{'GOOD'}{'R'}{'me'}};
			my $Dme_count = keys %{$$localised{'GOOD'}{'D'}{'me'}};
			my $Eme_count = keys %{$$localised{'GOOD'}{'E'}{'me'}};
			my $Qme_count = keys %{$$localised{'GOOD'}{'Q'}{'me'}};
			my $Hme_count = keys %{$$localised{'GOOD'}{'H'}{'me'}};		
			my $Kdi_count = keys %{$$localised{'GOOD'}{'K'}{'di'}};
			my $Rdi_count = keys %{$$localised{'GOOD'}{'R'}{'di'}};
			my $Mox_count = keys %{$$localised{'GOOD'}{'M'}{'ox'}};
			my $Ktr_count = keys %{$$localised{'GOOD'}{'K'}{'tr'}};
			my $Met_count = keys %{$$localised{'GOOD'}{'M'}{'me'}};		#va a contare quante modificazioni sono state assegnate in modo affidabile
			
			next if (($Methyl_all > 0) && (($Kme_count + $Rme_count + $Dme_count + $Eme_count + $Qme_count + $Hme_count) != $Methyl_all));
			next if (($Methyl > 0) && (($Kme_count + $Rme_count) != $Methyl));
			next if (($Dimethyl > 0) && (($Kdi_count + $Rdi_count) != $Dimethyl));
			next if (($Trimethyl > 0) && ($Ktr_count != $Trimethyl));				#ricerca errori=il numero di modifiche totali deve esse = somma delle singole modifiche

			next if (($Methyl4_all > 0) && (($Kme_count + $Rme_count + $Dme_count + $Eme_count + $Qme_count + $Hme_count) != $Methyl4_all));
			next if (($Methyl4 > 0) && (($Kme_count + $Rme_count) != $Methyl4));
			next if (($Dimethyl4 > 0) && (($Kdi_count + $Rdi_count) != $Dimethyl4));
			next if (($Trimethyl4 > 0) && ($Ktr_count != $Trimethyl4));
			next if ($Mox_count != $Oxidation);		
			next if (($Met4 > 0) && ($Met_count != $Met4));
			next if (($Met4ox > 0) && ($Mox_count != $Met4ox));			#se ci sono sbagli salta la riga e analizza quella dopo
										
			my ($best_msms_id) = ($1) if ($columns[$evidence_columns{$file}{"Best MS/MS"}] =~ /^(\d+)\;*.*/);
			my $delta_mass;
			my $methyl_diff_mass = 4.0221850755;

			my $no_m = $Modified_Sequence;
			$no_m =~ s/M//g;
			my $m_freq = length($Modified_Sequence)-length($no_m);
			
			# CALCULATE MASS SHIFT
			if ($Methyl_all > 0 || $Methyl > 0 || $Dimethyl > 0 || $Trimethyl > 0)
			{
			$delta_mass = ($Methyl_all * $methyl_diff_mass) + ($Methyl * $methyl_diff_mass)  + ($Dimethyl * $methyl_diff_mass * 2) + ($Trimethyl * $methyl_diff_mass * 3) + ($Oxidation * $methyl_diff_mass) + (($m_freq-$Oxidation) * $methyl_diff_mass);
			}
			elsif ($Methyl4_all > 0 || $Methyl4 > 0 || $Dimethyl4 > 0 || $Trimethyl4 > 0)
			{
			$delta_mass = 0-(($Methyl4_all * $methyl_diff_mass) + ($Methyl4 * $methyl_diff_mass)  + ($Dimethyl4 * $methyl_diff_mass * 2) + ($Trimethyl4 * $methyl_diff_mass * 3) + ($Met4 * $methyl_diff_mass) + ($Met4ox * $methyl_diff_mass));
			}
			
			foreach my $aa (keys %{$$localised{'GOOD'}})
			{
				foreach my $type (keys %{$$localised{'GOOD'}{$aa}})
				{
					foreach my $pos (keys %{$$localised{'GOOD'}{$aa}{$type}})
					{	
						my $line;
						
						next unless ($aa =~ /K|R|D|E|Q|H/);		#ignora gli aa diversi da questi
					
						my $modification_window = $$localised{'GOOD'}{$aa}{$type}{$pos}{'window'};
						my $residue_number		= $$localised{'GOOD'}{$aa}{$type}{$pos}{'residue'};
						my $state				= $$localised{'GOOD'}{$aa}{$type}{$pos}{'state'};

						$line .=  $state."\t";
						$line .=  $file."\t";
						$line .=  $columns[$evidence_columns{$file}{"id"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"Modified sequence"}]."\t";
						$line .=  $modification_window."\t";
						$line .=  $columns[$evidence_columns{$file}{"Missed cleavages"}]."\t";
						$line .=  $leading_proteins."\t";
						$line .=  $columns[$evidence_columns{$file}{"Gene Names"}]."\t";
						$line .=  "$aa\t$type\t$residue_number\t";
						$line .=  $columns[$evidence_columns{$file}{"PEP"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"Score"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"Delta score"}]."\t";
						$line .=  $experiment."\t";
						$line .=  $columns[$evidence_columns{$file}{"Raw file"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"m/z"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"Intensity"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"Retention time"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"MS/MS Scan Number"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"Best MS/MS"}]."\t";
						$line .=  $msms_info{$best_msms_id}{'Peak Coverage'}."\t";
						$line .=  $columns[$evidence_columns{$file}{"Sequence"}]."\t";
						$line .=  $columns[$evidence_columns{$file}{"Charge"}]."\t";
						$line .=  $msms_info{$best_msms_id}{'Peak List'}."\t";
						$line .=  $msms_info{$best_msms_id}{'Peak List Corrected'}."\t";
						$line .=  $delta_mass."\t";
						$line .=  $columns[$evidence_columns{$file}{"Protein Names"}];
							
						print OUT_REDUNDANT "$line\n" unless ($line =~ /^$/);
					}
				}
			}
		}
	}
	close EVIDENCE;
	close OUT_REDUNDANT;
}

###################### DOUBLETS SEARCH ###############################################
my (%results_columns,%reference,%raw_file_limits,%scan_results);
my (%results_columns);
my $line;

foreach my $file (@files)	#READ INPUT FILES (msmsScans;combined_all_red)
{
my $subfolder = $file."_".$date."_".$hour.$min;
open IN, $output_folder.$subfolder.'/peptides_summary.txt' or die "COULD NOT FIND RESULTS";
while ($line=<IN>)
{
	chomp($line);
	my (@columns) = split(/\t/,$line);
	
	#DATI PRESI DA ALL REDUNDANT 
	if ($line =~ /^State/)
	{
		(%results_columns) = %{&maxquant::columns(\@columns)};
	}
	else
	{
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'dM'} = $columns[$results_columns{"Delta Mass"}];
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'State'} = $columns[$results_columns{"State"}];
		#$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Precursor Intensity'} = $columns[$results_columns{"Intensity"}];
				
		#$raw_file_limits{$columns[$results_columns{"Raw file"}]}{"Min"} = $columns[$results_columns{"Retention time"}] if (!defined($raw_file_limits{$columns[$results_columns{"Raw file"}]}{"Min"}) || (defined($raw_file_limits{$columns[$results_columns{"Raw file"}]}{"Min"}) && $raw_file_limits{$columns[$results_columns{"Raw file"}]}{"Min"} > $columns[$results_columns{"Retention time"}]));
		#$raw_file_limits{$columns[$results_columns{"Raw file"}]}{"Max"}	= $columns[$results_columns{"Retention time"}] if (!defined($raw_file_limits{$columns[$results_columns{"Raw file"}]}{"Max"}) || (defined($raw_file_limits{$columns[$results_columns{"Raw file"}]}{"Max"}) && $raw_file_limits{$columns[$results_columns{"Raw file"}]}{"Max"} < $columns[$results_columns{"Retention time"}]));	
	
$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Sequence'} 				= $columns[$results_columns{"Sequence"}];
$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Modified sequence'}		= $columns[$results_columns{"Modified sequence"}];
$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Modification window'}		= $columns[$results_columns{"Modification Window"}];
	
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Residue'} 		= $columns[$results_columns{"Residue"}];
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Modification'} 	= $columns[$results_columns{"Modification"}];
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Position'} 		= $columns[$results_columns{"Position in protein"}];
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Proteins'} 		= $columns[$results_columns{"Leading Proteins"}];
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Score'} 		= $columns[$results_columns{"Score"}];
		my @gene_names = split(";", $columns[$results_columns{"Gene Names"}]);
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Gene'} 		= $gene_names[0];
		$reference{$columns[$results_columns{"Dataset"}]}{$columns[$results_columns{"Raw file"}]}{$columns[$results_columns{"MS/MS Scan Number"}]}{'Protein name'} 		= $columns[$results_columns{"Protein Names"}];
	
	
	}
}
close IN;

print "stampo ".$file."_msmsScans\n";

	my (%msmsScan_columns,$i);

	#DATI PRESI DA MSMSSCANS 
	open MSMSScans, $dir.$file."/msmsScans.txt" or die "Could not open ".$file;
	open OUT, '>'.$output_folder.$subfolder.'/msms_scans.txt';
	while (<MSMSScans>)
	{
		chomp($_);
		my (@columns) = split(/\t/,$_);

		if ($_ =~ /^Raw/)
		{
			(%msmsScan_columns) = %{&maxquant::columns(\@columns)};
			print OUT "$_\n";
		}
		else
		{
			my $scan_number = $columns[$msmsScan_columns{"Scan number"}];
			my $rawfile     = $columns[$msmsScan_columns{"Raw file"}];
			next if ($columns[$msmsScan_columns{"Charge"}] < 2);
			#next unless (defined($raw_file_limits{$columns[$msmsScan_columns{"Raw file"}]}{"Min"}) && defined($raw_file_limits{$columns[$msmsScan_columns{"Raw file"}]}{"Max"}) && ($columns[$msmsScan_columns{"Retention time"}] <= $raw_file_limits{$columns[$msmsScan_columns{"Raw file"}]}{"Max"}) && ($columns[$msmsScan_columns{"Retention time"}] >= $raw_file_limits{$columns[$msmsScan_columns{"Raw file"}]}{"Min"}));
			
			$i++;

			$scan_results{$rawfile}{'details'}{$scan_number}{'Charge'}					= $columns[$msmsScan_columns{"Charge"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'Identified'} 				= $columns[$msmsScan_columns{"Identified"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'m/z'} 					= $columns[$msmsScan_columns{"m/z"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'Raw file'} 				= $columns[$msmsScan_columns{"Raw file"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'Retention time'} 			= $columns[$msmsScan_columns{"Retention time"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'Precursor intensity'}		= $columns[$msmsScan_columns{"Precursor intensity"}];
			###parte aggiunta da me###
			$scan_results{$rawfile}{'details'}{$scan_number}{'Sequence'}				= $columns[$msmsScan_columns{"Sequence"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'Modified sequence'}		= $columns[$msmsScan_columns{"Modified sequence"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'Score'}					= $columns[$msmsScan_columns{"Score"}];	
			$scan_results{$rawfile}{'details'}{$scan_number}{'Modifications'}			= $columns[$msmsScan_columns{"Modifications"}];	
			$scan_results{$rawfile}{'details'}{$scan_number}{'Proteins'}				= $columns[$msmsScan_columns{"Proteins"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'Gene'}					= $columns[$msmsScan_columns{"Gene Names"}];
			$scan_results{$rawfile}{'details'}{$scan_number}{'Protein name'}			= $columns[$msmsScan_columns{"Protein Names"}];
			
			$scan_results{$rawfile}{'details'}{$scan_number}{'Index'} 	= $i;
			$scan_results{$rawfile}{'index'}{$i} 	= $scan_number;

			print OUT "$_\n"; #if (defined($raw_file_limits{$columns[$msmsScan_columns{"Raw file"}]}{"Min"}) && defined($raw_file_limits{$columns[$msmsScan_columns{"Raw file"}]}{"Max"}) && ($columns[$msmsScan_columns{"Retention time"}] <= $raw_file_limits{$columns[$msmsScan_columns{"Raw file"}]}{"Max"}) && ($columns[$msmsScan_columns{"Retention time"}] >= $raw_file_limits{$columns[$msmsScan_columns{"Raw file"}]}{"Min"}));
		}
	}
	close MSMSScans;
	close OUT;

	
print "stampo ".$file."_doublets_ALL\n";

my %pep_found;
my %not_found;
my $doublet;
my $counter = 1;
##prepare doublets_all
open OUT, '>'.$output_folder.$subfolder.'/doublets_ALL.txt';
my $header = "Raw file\tH Scan number\tL Scan number\tH Mass/charge\tL Mass/charge\tDelta Mass (observed)\tDelta Mass (theoretic)\tH RT\tL RT\tH Charge\tL Charge\tH ID\tL ID\tRaw_Scan H\tRaw_Scan L\tH Intensity\tL Intensity\tScore H\tScore L\t";
$header .= "Sequence H\tSequence L\tModified Sequence H\tModified Sequence L\tWindow H\tWindow L\tMismatch\tResidue H\tResidue L\tPosition H\tPosition L\tModification H\tModification L\t";
$header .= "HMet\tHMetOx\tMet\tMetOx\tProtein H\tProtein L\tGene H\tGene L\tProtein Name H\tProtein Name L\tID number";
print OUT "$header\n";

foreach my $rawfile (keys %{$reference{$file}})
	{
	foreach my $scan_number (keys %{$reference{$file}{$rawfile}})
		{
			my $delta_mass		= $reference{$file}{$rawfile}{$scan_number}{'dM'};
			my $state			= $reference{$file}{$rawfile}{$scan_number}{'State'};	
			
			#extracts data from msmsScans
			my $index 			= $scan_results{$rawfile}{'details'}{$scan_number}{'Index'};
			my $charge 			= $scan_results{$rawfile}{'details'}{$scan_number}{'Charge'};
				next if ($charge eq "");	#filtro
			
			my $mz 				= $scan_results{$rawfile}{'details'}{$scan_number}{'m/z'};
			my $Identified 		= $scan_results{$rawfile}{'details'}{$scan_number}{'Identified'};
			my $RT 				= $scan_results{$rawfile}{'details'}{$scan_number}{'Retention time'};
			my $intensity		= $scan_results{$rawfile}{'details'}{$scan_number}{'Precursor intensity'};

			my $mod_seq  =  $scan_results{$rawfile}{'details'}{$scan_number}{'Modified sequence'}; #tries to get data from msmsScans

			my $pep_seq;
			my @res;
			my @mod;
			my @pos;
			my @windows;
			my $prot;
			my $gene		= $reference{$file}{$rawfile}{$scan_number}{'Gene'};
			my $prot_name	= $reference{$file}{$rawfile}{$scan_number}{'Protein name'};
			my $score;
			my $met; 
			my $metox;
			
			if ($mod_seq eq " " || $mod_seq eq "" || !defined $mod_seq) #data from allred if for some reason msmsScans does not contain the data ($mod_seq eq " " || $mod_seq eq "" || !defined $mod_seq)
			{
				$pep_seq			= $reference{$file}{$rawfile}{$scan_number}{'Sequence'};	
				$mod_seq			= $reference{$file}{$rawfile}{$scan_number}{'Modified sequence'};
				$res[0]				= $reference{$file}{$rawfile}{$scan_number}{'Residue'};
				$mod[0] 			= $reference{$file}{$rawfile}{$scan_number}{'Modification'} ;
				$pos[0]				= $reference{$file}{$rawfile}{$scan_number}{'Position'};
				$prot				= $reference{$file}{$rawfile}{$scan_number}{'Proteins'};
				$score				= $reference{$file}{$rawfile}{$scan_number}{'Score'};
				$windows[0]			= $reference{$file}{$rawfile}{$scan_number}{'Modification window'};
				
				$mod_seq =~ s/\_//g;
				#$mod_seq =~ s/\(ac\)//g;
			}
			else         #data from msmsscans
			{
				$pep_seq  	= 	$scan_results{$rawfile}{'details'}{$scan_number}{'Sequence'};  
				$score 		=  	$scan_results{$rawfile}{'details'}{$scan_number}{'Score'};
				$prot		= 	$scan_results{$rawfile}{'details'}{$scan_number}{'Proteins'};
				
				#mapping modification (in msmsScans there is no mapping)
				my @leading_proteins = split(';',$prot);
				my @position_on_protein;	
				$mod_seq =~ s/\_//g;
				my $position_on_pep = 0;
				my $n = 0;
				my $offset = 0;
				while ($position_on_pep != -1)
					{
					$position_on_pep = index($mod_seq,'(',$position_on_pep+1);
					if ( substr($mod_seq,$position_on_pep+1,2) =~ /me|di|tr/ && substr($mod_seq,$position_on_pep-1,1) ne 'M' && $position_on_pep != -1)
						{
						$res[$n] = substr($mod_seq,$position_on_pep-1,1);
						$mod[$n] = substr($mod_seq,$position_on_pep+1,2);
						$pos[$n] = $position_on_pep - $offset;
						$n++;
						}
					$offset += 4;
					}
				
				my $pep_seq_cleared = $pep_seq;
				$pep_seq_cleared =~ s/\I|L/X/g;
				$position_on_protein[0] = index($fasta_database{$leading_proteins[0]}{'IL_cleaned'},$pep_seq_cleared);
				#determine sequence windows
				for ($n=0; $n<scalar(@pos);$n++)
					{
					$pos[$n] += $position_on_protein[0];
					#my $buffer = 'zzzzzzzz';
					#$windows[$n] = substr($buffer.$fasta_database{$leading_proteins[0]}{'Raw'}.$buffer, $pos[$n], 15);
					
					my $buffer = 'z' x $Mod_Window;
					my $sequence = $buffer.$fasta_database{$leading_proteins[0]}{'Raw'}.$buffer;
					my $window_start   = ($Mod_Window+$pos[$n]-(($Mod_Window-1)/2));	
					$windows[$n]	   = substr($sequence,$window_start-1,$Mod_Window);
					}
			}
			
			#controllo sulle metionine
			#next if ($state =~ /H/ && $mod_seq =~ /M/ && $mod_seq !~ /M\(me\)/ && $mod_seq !~ /M\(ox\)/);
			#next if ($state =~ /L/ && $mod_seq =~ /M\(me\)/);
			#presence of methynine (pesanti/leggere/ossidate/ridotte)
			if ($mod_seq =~ /M/)
				{$met = '+';}
			elsif ($mod_seq !~ /^$/)
				{$met = '-';}
				
			if ($mod_seq =~ /ox/)
				{$metox = '+';}
			elsif ($mod_seq !~ /^$/)
				{$metox = '-';}
			
	####counter check - this will tell the script if the object has a counterpart or not
	 my $counter_check = $counter;
			
			for (my $j=($index-999);$j<($index+999);$j++)
			#foreach my $j (keys %{$scan_results{$rawfile}{'index'}})
				{
				my $new_scan = $scan_results{$rawfile}{'index'}{$j};
				my $new_state		= $reference{$file}{$rawfile}{$new_scan}{'State'};
				next if ($new_state eq $state);
					#extracts data from msmsScans
				my $new_charge 		= $scan_results{$rawfile}{'details'}{$new_scan}{'Charge'};
				my $new_mz 			= $scan_results{$rawfile}{'details'}{$new_scan}{'m/z'};
				my $new_Identified 	= $scan_results{$rawfile}{'details'}{$new_scan}{'Identified'};
				my $new_RT 			= $scan_results{$rawfile}{'details'}{$new_scan}{'Retention time'};
				my $new_intensity	= $scan_results{$rawfile}{'details'}{$new_scan}{'Precursor intensity'};
				my $new_mod_seq  	= $scan_results{$rawfile}{'details'}{$new_scan}{'Modified sequence'}; #tries to get data from msmsScans

				my $new_pep_seq;
				my @new_res;
				my @new_mod;
				my @new_pos;
				my @new_windows;
				my $new_prot;
				my $new_gene		= $reference{$file}{$rawfile}{$new_scan}{'Gene'};
				my $new_prot_name	= $reference{$file}{$rawfile}{$new_scan}{'Protein name'};
				my $new_score;
				my $new_met;
				my $new_metox;
			
				if ($new_mod_seq eq " " || $new_mod_seq eq "" || !defined $new_mod_seq) #data from allred if or some reason msmsScans does not contain the data
					{
					$new_pep_seq			= $reference{$file}{$rawfile}{$new_scan}{'Sequence'}; 
				
					$new_mod_seq			= $reference{$file}{$rawfile}{$new_scan}{'Modified sequence'};
					$new_res[0]				= $reference{$file}{$rawfile}{$new_scan}{'Residue'};
					$new_mod[0]				= $reference{$file}{$rawfile}{$new_scan}{'Modification'};
					$new_pos[0]				= $reference{$file}{$rawfile}{$new_scan}{'Position'};
					$new_prot				= $reference{$file}{$rawfile}{$new_scan}{'Proteins'};
					$new_score				= $reference{$file}{$rawfile}{$new_scan}{'Score'};
					$new_windows[0]			= $reference{$file}{$rawfile}{$new_scan}{'Modification window'};
					
					$new_mod_seq =~ s/\_//g;
					#$new_mod_seq =~ s/\(ac\)//g;
					}
				else  #data from msmsscans
					{
					$new_pep_seq  	=  $scan_results{$rawfile}{'details'}{$new_scan}{'Sequence'}; 	
					$new_score 		=  $scan_results{$rawfile}{'details'}{$new_scan}{'Score'};
					$new_prot		=  $scan_results{$rawfile}{'details'}{$new_scan}{'Proteins'};	
				
				#mapping modification (in msmsScans there is no mapping)
					my @leading_proteins = split(';',$new_prot);
					my @position_on_protein;	
					$new_mod_seq =~ s/\_//g;
					my $position_on_pep = 0;
					my $n = 0;
					my $offset = 0;
					while ($position_on_pep != -1)
						{
						$position_on_pep = index($new_mod_seq,'(',$position_on_pep+1);
						if ( substr($new_mod_seq,$position_on_pep+1,2) =~ /me|di|tr/ && substr($new_mod_seq,$position_on_pep-1,1) ne 'M' && $position_on_pep != -1)
							{
							$new_res[$n] = substr($new_mod_seq,$position_on_pep-1,1);
							$new_mod[$n] = substr($new_mod_seq,$position_on_pep+1,2);
							$new_pos[$n] = $position_on_pep - $offset;
							$n++;
							}
						$offset += 4;
						}
					my $pep_seq_cleared = $new_pep_seq;
					$pep_seq_cleared =~ s/\I|L/X/g;
					my $ref_seq = $fasta_database{$leading_proteins[0]}{'IL_cleaned'};
					$position_on_protein[0] = index($ref_seq,$pep_seq_cleared);
				#determine seq windows
					for ($n=0; $n<scalar(@new_pos);$n++)
						{	
						$new_pos[$n] += $position_on_protein[0];
						#my $buffer = 'zzzzzzzz';
						#$new_windows[$n] = substr($buffer.$fasta_database{$leading_proteins[0]}{'Raw'}.$buffer, $new_pos[$n], 15);
						
						my $buffer = 'z' x $Mod_Window;
						my $sequence = $buffer.$fasta_database{$leading_proteins[0]}{'Raw'}.$buffer;
						my $window_start   = ($Mod_Window+$new_pos[$n]-(($Mod_Window-1)/2));	
						$new_windows[$n]	   = substr($sequence,$window_start-1,$Mod_Window);
						}
						
					}
					
			#controllo metionine
			#	next if ($new_state =~ /H/ && $new_mod_seq =~ /M/ && $new_mod_seq !~ /M\(me\)/ && $new_mod_seq !~ /M\(ox\)/);
			#	next if ($new_state =~ /L/ && $new_mod_seq =~ /M\(me\)/);	
			#presence of methyonines
				if ($new_mod_seq =~ /M/)
					{$new_met = '+';}
				elsif ($new_mod_seq !~ /^$/)
					{$new_met = '-';}
				
				if ($new_mod_seq =~ /ox/)
					{$new_metox = '+';}
				elsif ($new_mod_seq !~ /^$/)
					{$new_metox = '-';}
				
			#modified sequence mismatch
				$mod_seq =~ s/M\(me\)/M/g;
				$mod_seq =~ s/M\(ox\)/M/g;
				$new_mod_seq =~ s/M\(me\)/M/g;
				$new_mod_seq =~ s/M\(ox\)/M/g;
				my $mis;
				$mis = "+" if ( $mod_seq ne $new_mod_seq && $mod_seq !~ /^$/ && $new_mod_seq !~ /^$/ );
				
				if ( ($state =~ /H/) && ($mz > $new_mz) && (abs($RT-$new_RT) < $RT_threshold) && ($charge == $new_charge) && (abs((($mz-$new_mz)*$charge)+$delta_mass) < $mass_threshold))
					{
					my $raw_scan_H = $rawfile."_".$scan_number;
					my $raw_scan_L = $rawfile."_".$new_scan;
					foreach my $key (keys @pos)
						{
						$doublet = "$rawfile\t$scan_number\t$new_scan\t$mz\t$new_mz\t".(abs($new_mz-$mz)*$charge)."\t".abs($delta_mass)."\t$RT\t$new_RT\t$charge\t$new_charge\t$Identified\t$new_Identified\t$raw_scan_H\t$raw_scan_L\t$intensity\t$new_intensity\t$score\t$new_score\t";
						$doublet .= "$pep_seq\t$new_pep_seq\t$mod_seq\t$new_mod_seq\t$windows[$key]\t$new_windows[$key]\t$mis\t$res[$key]\t$new_res[$key]\t$pos[$key]\t$new_pos[$key]\t$mod[$key]\t$new_mod[$key]\t$met\t$metox\t$new_met\t$new_metox\t$prot\t$new_prot\t$gene\t$new_gene\t$prot_name\t$new_prot_name\t$counter";
						print OUT "$doublet\n";
						$counter++;
						}
					}
				elsif (($state =~ /L/) && ($mz < $new_mz) && (abs($RT-$new_RT) < $RT_threshold) && ($charge == $new_charge) && (abs((($mz-$new_mz)*$charge)+$delta_mass) < $mass_threshold))
					{
					my $raw_scan_H = $rawfile."_".$new_scan;
					my $raw_scan_L = $rawfile."_".$scan_number;
					foreach my $key (keys @pos)
						{						
						$doublet = "$rawfile\t$new_scan\t$scan_number\t$new_mz\t$mz\t".(abs($new_mz-$mz)*$charge)."\t".abs($delta_mass)."\t$new_RT\t$RT\t$new_charge\t$charge\t$new_Identified\t$Identified\t$raw_scan_H\t$raw_scan_L\t$new_intensity\t$intensity\t$new_score\t$score\t";
						$doublet .= "$new_pep_seq\t$pep_seq\t$new_mod_seq\t$mod_seq\t$new_windows[$key]\t$windows[$key]\t$mis\t$new_res[$key]\t$res[$key]\t$new_pos[$key]\t$pos[$key]\t$new_mod[$key]\t$mod[$key]\t$new_met\t$new_metox\t$met\t$metox\t$new_prot\t$prot\t$new_gene\t$gene\t$new_prot_name\t$prot_name\t$counter";
						print OUT "$doublet\n";
						$counter++;
						}
					}
			
				}
				
			my $sites = "";
			if ($counter_check == $counter)		#if the script doesn't find doublets, the counter in not incremented @ lines 675/687
				{
				foreach my $key (keys @pos)
					{
					$sites .= $res[$key]."_".$pos[$key]."_".$mod[$key].";";
					}
				my $object = "$state\t$Identified\t$rawfile\t$scan_number\t$mz\t$RT\t$charge\t".($delta_mass/$charge)."\t$intensity\t$score\t$pep_seq\t$mod_seq\t$sites\t$gene\t$prot_name";
				push (@{$not_found{$mod_seq}}, $object);
				}
			else
				{
				$pep_found{$mod_seq} = 1;
				}
		}
	
	}

close OUT;

############### PEPETIDES NOT FOUND #############
print "stampo ".$file."_NotFound\n";

my $novel;
##prepare notfound
open NOTFOUND, '>'.$output_folder.$subfolder.'/NotFound.txt';
my $header = "NOVEL\tState\tIdentified\tRawfile\tScan Number\tMass\tRT\tCharge\tDelta m/z (theoretic)\tIntensity\tScore\tSequence\tModified Sequence\tModifications Summary\tGene\tProtein Name";
print NOTFOUND "$header\n";

foreach my $mod_seq (keys %not_found)
	{
	foreach my $object (@{$not_found{$mod_seq}}) 
		{
		$novel = "";
		if (!exists $pep_found{$mod_seq})	#check novelty
			{
			$novel = "+";
			$pep_found{$mod_seq} = 1;
			}
		
		print NOTFOUND "$novel\t$object\n";
		
		}
	
	}
close NOTFOUND;

######################################
####### NON REDUNDANT SILAC ##########
print "stampo ".$file."_doublets_NR\n";

my $line;
my $header = "Raw file\tH Scan number\tL Scan number\tH Mass/charge\tL Mass/charge\tDelta Mass (observed)\tDelta Mass (theoretic)\tH RT\tL RT\tH Charge\tL Charge\tH ID\tL ID\tRaw_Scan H\tRaw_Scan L\tH Intensity\tL Intensity\tScore H\tScore L\t";
$header .= "Sequence H\tSequence L\tModified Sequence H\tModified Sequence L\tWindow H\tWindow L\tMismatch\tResidue H\tResidue L\tPosition H\tPosition L\tModification H\tModification L\t";
$header .= "HMet\tHMetOx\tMet\tMetOx\tProtein H\tProtein L\tGene H\tGene L\tProtein Name H\tProtein Name L\tID number";
my $rawfile;
my %doublets_columns;
my %nr_lines;
my @momentaneo;
my @columns;

#input
open REDUNDANT, $output_folder.$subfolder.'/doublets_ALL.txt';
while ($line = <REDUNDANT> )
	{
	chomp ($line);
	@columns = split(/\t/,$line);
	
	if ($line =~ /^Raw/)
		{
		(%doublets_columns) = %{&maxquant::columns(\@columns)};
		}
	else
		{
		next if $columns[$doublets_columns{"Mismatch"}] eq "+";		#FILTRO
		
		my $seq;
		my $res;
		my $pos;
		my $mod;
		
		if ($columns[$doublets_columns{"Modified Sequence H"}] !~ /^$/)	#extracts data
			{
			$seq	=	$columns[$doublets_columns{"Modified Sequence H"}];
			$res	=	$columns[$doublets_columns{"Residue H"}];
			$pos	=	$columns[$doublets_columns{"Position H"}];
			$mod	=	$columns[$doublets_columns{"Modification H"}];
			}
		else
			{
			$seq 	=	$columns[$doublets_columns{"Modified Sequence L"}];
			$res	=	$columns[$doublets_columns{"Residue L"}];
			$pos	=	$columns[$doublets_columns{"Position L"}];
			$mod	=	$columns[$doublets_columns{"Modification L"}];
			}
		
		if (exists $nr_lines{$seq}{$res}{$pos}{$mod})	#checks if a doublet already exists
			{
			@momentaneo = split ('\t', $nr_lines{$seq}{$res}{$pos}{$mod});
			#priority
			my $ID = 0;
			$ID++ if ($momentaneo[$doublets_columns{"L ID"}] =~ /\+/);
			$ID += 2 if ($momentaneo[$doublets_columns{"H ID"}] =~ /\+/);
			my $new_ID = 0;
			$new_ID++ if ($columns[$doublets_columns{"L ID"}] =~ /\+/);
			$new_ID += 2 if ($columns[$doublets_columns{"H ID"}] =~ /\+/);
			#scores
			my $average_score =			$momentaneo[$doublets_columns{"Score H"}]+$momentaneo[$doublets_columns{"Score L"}];
			my $new_average_score = 	$columns[$doublets_columns{"Score H"}]+$columns[$doublets_columns{"Score L"}];
			
			if ($new_ID>$ID || ($ID==$new_ID && $average_score<$new_average_score))		#compares priority and score
				{
				$nr_lines{$seq}{$res}{$pos}{$mod} = $line.';'.$momentaneo[$doublets_columns{"ID number"}];	 #adds ID of the discarded dobulet
				}
			else
				{
				$nr_lines{$seq}{$res}{$pos}{$mod} .= ';'.$columns[$doublets_columns{"ID number"}];	 #adds ID of the discarded dobulet
				}
			}
		else 
			{
			$nr_lines{$seq}{$res}{$pos}{$mod} = $line;
			}
		}
	}


	
close REDUNDANT;

#print output
open NR, '>'.$output_folder.$subfolder.'/doublets_NR.txt';
print NR "$header\n";
foreach my $seq (keys %nr_lines)
	{
		foreach my $res (sort (keys %{$nr_lines{$seq}}))
			{
				foreach my $pos (keys %{$nr_lines{$seq}{$res}})
					{
						foreach my $mod (keys %{$nr_lines{$seq}{$res}{$pos}})
							{
								print NR "$nr_lines{$seq}{$res}{$pos}{$mod}\n";
							}
					}
			}
	}
	
close NR;

############PRINTING PARAMETERS#########
open par, '>'.$output_folder.$subfolder."/parameters.txt";
print par "Run date (yy_mm_dd):\t$date\nRetention Time threshold (min):\t$RT_threshold\nMass Error threshold (Da):\t$mass_threshold\nMinimum Score:\t$Score_Min\n";
print par "Minimum Delta Score:\t$Delta_Min\nMinimum Localization Probability:\t$LcPrb_Min\nModification Window:\t$Mod_Window\nFASTA Database:\t$database\n";
print par "Input Datasets:";
foreach my $file (@files)
	{ print par "\t$file";	}


}
print "\nDONE\n";