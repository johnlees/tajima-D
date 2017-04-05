/*
 * File: tajima.cpp
 *
 * Main control loop for tajima
 *
 */

#include "tajima.hpp"

// Constants
const std::string VERSION = "0.1";
//    Default options

int main (int argc, char *argv[])
{
   // Read cmd line options and process
   // Read in all variant lines, store positions of 1s in vector for each line
   // Create N Sample objects, with var size
   // For each variant line, add any ones into relevant sample at correct
   // positon
   // Using full sample objects, calculate Tajima's D
   //
   // example command to generate input
   // bcftools view -f PASS -r FM211187:186-1547 /lustre/scratch108/bacteria/jl11/mappings/23FSpn_new/recalibrated.snps.indels.vcf.gz | bcftools norm -m - | bcftools query -f '[%GT,]\n' | sed 's/,$//'

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: tajima --snps snps.csv\n\n"
         << "For full option details run tajima -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   // Program description
   if (vm.count("verbose"))
   {
      std::cerr << "tajima: calculate tajima's D from bcftools query output\n";
   }

   // Read in all variant lines
   if (vm.count("verbose"))
   {
      std::cerr << "Reading all variants in..." << std::endl;
   }
   std::vector<std::vector<int>> variants;
   if (vm.count("snps"))
   {
      std::ifstream file_in;
      file_in.open(vm["snps"].as<std::string>().c_str());
      while (file_in)
      {
         variants.push_back(readCsvLine(file_in));
      }
   }
   else
   {
      throw std::logic_error("--snps option is compulsory");
   }

   // Set up vector of samples
   size_t num_samples = variants[0].size();
   std::vector<Sample> samples;
   samples.reserve(num_samples);
   for (size_t i = 0; i < num_samples; i++)
   {
      Sample s(variants.size());
      samples.push_back(s);
   }

   // Transpose variants into samples
   if (vm.count("verbose"))
   {
      std::cerr << "Transposing..." << std::endl;
   }
   for (size_t var_it = 0; var_it < variants.size(); ++var_it)
   {
      for (size_t samp_it = 0; samp_it < variants[var_it].size(); ++samp_it)
      {
         if (variants[var_it][samp_it] == 1)
         {
            samples[samp_it].add_to_seq(1, var_it);
         }
      }
   }

   // Calculate D
   // https://en.wikipedia.org/wiki/Tajima's_D#Mathematical_details
   if (vm.count("verbose"))
   {
      std::cerr << "Calculating D..." << std::endl;
   }
   double d_sum = 0;
   int d_tot = 0;
   for (size_t i = 0; i < samples.size(); i++)
   {
      for (size_t j = i+1; j < samples.size(); j++)
      {
         d_sum += dot(samples[i].full_seq() - samples[j].full_seq(), samples[i].full_seq() - samples[j].full_seq());
         d_tot++;
      }
   }

   double k_hat = d_sum/d_tot;
   double S = variants.size();

   // Override segregating sites. If multiallelics are present the default will
   // overcount them
   if (vm.count("seg_sites"))
   {
      S = vm["seg_sites"].as<double>();
   }

   if (vm.count("verbose"))
   {
      std::cerr << "k_hat: " << k_hat << std::endl;
      std::cerr << "S: " << S << std::endl;
   }

   double a1 = 0, a2 = 0;
   for (size_t i = 1; i <= num_samples - 1; i++)
   {
      a1 += pow(i,-1);
      a2 += 1/pow(i,2);
   }

   double b1 = ((double)num_samples + 1)/(3*(num_samples - 1));
   double b2 = 2*(pow(num_samples,2) + num_samples + 3)/(9*num_samples*(num_samples - 1));

   double c1 = b1 - 1/a1;
   double c2 = b2 - (num_samples + 2)/(a1*num_samples) + a2/pow(a1,2);

   double e1 = c1/a1;
   double e2 = c2/(pow(a1, 2) + a2);

   double se = pow(e1*S + e2*S*(S-1), 0.5);

   double D = (k_hat - S/a1)/se;

   std::cout << D << std::endl;

   if (vm.count("verbose"))
   {
      std::cerr << "Done.\n";
   }
}

std::vector<int> readCsvLine(std::istream& is)
{
   std::vector<int> variant;
   std::string line;
   std::getline(is, line);

   std::stringstream line_stream(line);
   std::string value;

   while(std::getline(line_stream, value, ','))
   {
      if (value == "1")
      {
         variant.push_back(1);
      }
      else
      {
         variant.push_back(0);
      }
   }
   return variant;
}

