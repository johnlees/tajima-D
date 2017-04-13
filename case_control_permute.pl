#!/usr/bin/perl -w

use strict;
use warnings;

use File::Temp qw/ :POSIX /;
use List::Util 'shuffle';

sub run_tajima
{
   my ($vcf_file, $region, $samples, $exclude) = @_;

   my $norm_tmp = tmpnam();
   my $norm_err_tmp = tmpnam();
   my $input_tmp = tmpnam();

   my $bcftools_command = "bcftools view -f PASS -r $region";
   if ($exclude)
   {
      $bcftools_command .= " -s^" . join(",", @$samples);
   }
   else
   {
      $bcftools_command .= " -s" . join(",", @$samples);
   }

   $bcftools_command .= " -c 1 $vcf_file | bcftools norm -m - > $norm_tmp 2> $norm_err_tmp";
   system($bcftools_command);

   open(ERR, $norm_err_tmp) || die("Could not open $norm_err_tmp\n");
   my $line_in = <ERR>;
   chomp $line_in;
   $line_in =~ m/^Lines total\/modified\/skipped:\s+(\d+)\/(\d+)\/(\d+)$/;

   my $seg_sites = $1 - $3;
   my $tajima_D = 0;
   if ($seg_sites > 0)
   {
      if ($seg_sites < 3)
      {
         print STDERR "low seg sites: $seg_sites\n";
      }

      my $bcftools_command2 = "bcftools view -e 'ALT[*]==\"*\"' $norm_tmp | bcftools query -f '[%GT,]\n' | sed 's/,\$//' > $input_tmp";
      system($bcftools_command2);

      $tajima_D = `~/installations/tajima-D/tajima -s $input_tmp --seg_sites $seg_sites`;
      chomp $tajima_D;
   }

   unlink $norm_tmp, $norm_err_tmp, $input_tmp;
   return($tajima_D)
}

my $vcf_file = $ARGV[0];
my $region = $ARGV[1];
my $pheno_file = $ARGV[2];
my $permutations = $ARGV[3];

my $num_cases = 0;
my $num_controls = 0;

# Read in cases and controls
open(PHENO, $pheno_file) || die("Could not open $pheno_file\n");

my @control_samples;
my @case_samples;
my @all_samples;

while (my $line_in = <PHENO>)
{
   chomp $line_in;
   my ($fid, $iid, $pheno) = split("\t", $line_in);

   if ($pheno == 0)
   {
      push(@control_samples, $iid);
      $num_controls++;
   }
   else
   {
      push(@case_samples, $iid);
      $num_cases++;
   }

   push(@all_samples, $iid);
}

my @mixed_samples = shuffle(@all_samples);

my $case_D = run_tajima($vcf_file, $region, \@case_samples, 0);
my $control_D = run_tajima($vcf_file, $region, \@control_samples, 0);
my $diff_D = $case_D - $control_D;

my @D_permute;
for (my $i = 0; $i < $permutations; $i++)
{
   # Reservoir sampler
   my @case_array;
   for (my $j = 0; $j < $num_cases; $j++)
   {
      push(@case_array, $mixed_samples[$j]);
   }

   for (my $j = $num_cases; $j < $num_cases + $num_controls; $j++)
   {
      my $idx = int(rand($j));
      if ($idx < $num_cases)
      {
         $case_array[$idx] = $mixed_samples[$j]
      }
   }

   my $permuted_diff = run_tajima($vcf_file, $region, \@case_array, 0) - run_tajima($vcf_file, $region, \@case_array, 1);
   push(@D_permute, $permuted_diff);
}

my @D_permute_order = sort { $a <=> $b } @D_permute;
my $idx = 0;
foreach my $permuted_D (@D_permute_order)
{
   if ($diff_D < $permuted_D)
   {
      last;
   }
   $idx++;
}

if ($idx > scalar(@D_permute_order)/2)
{
   $idx = scalar(@D_permute_order) - $idx;
}

my $pval = 2*$idx/scalar(@D_permute_order);

print join("\t", $case_D, $control_D, $diff_D, $pval) . "\n";

exit(0);

