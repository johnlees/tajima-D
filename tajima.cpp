/*
 * File: tajima.cpp
 *
 * Main control loop for tajima
 *
 */

#include "tajima.hpp"

// Constants
const std::string VERSION = "0.3";
//    Default options

int main (int argc, char *argv[])
{
   // example command to generate input
   // bcftools view -f PASS -r <region> <vcf_file> | bcftools norm -m - | bcftools view -e 'ALT[*]=="*"' | bcftools query -f '%POS,[%GT,]\n' | sed 's/,$//' > input.csv

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
   std::vector<long int> positions;
   if (vm.count("snps"))
   {
      std::ifstream file_in;
      file_in.open(vm["snps"].as<std::string>().c_str());
      while (file_in)
      {
         std::string line;
         std::getline(file_in, line);
         if (file_in.good())
         {
            std::tuple<long int,std::vector<int>> variant = readCsvLine(line);
            positions.push_back(std::get<0>(variant));
            variants.push_back(std::get<1>(variant));
         }
         else
         {
            break;
         }
      }
   }
   else
   {
      throw std::logic_error("--snps option is compulsory");
   }

   // Transpose variants into samples
   if (vm.count("verbose"))
   {
      std::cerr << "Transposing..." << std::endl;
   }

   arma::mat alignment(variants[0].size(), variants.size(), arma::fill::zeros);
   for (size_t var_it = 0; var_it < variants.size(); ++var_it)
   {
      for (size_t samp_it = 0; samp_it < variants[var_it].size(); ++samp_it)
      {
         if (variants[var_it][samp_it] == 1)
         {
            alignment(samp_it, var_it) = 1;
         }
      }
   }

   if (vm.count("verbose"))
   {
      std::cerr << "Calculating D..." << std::endl;
   }

   size_t num_samples = alignment.n_rows;
   arma::mat pairwise_distances = dist_mat(alignment);
   // Run permutations subsampling from collection to generate null D
   if (vm.count("null_subsamples") && vm.count("null"))
   {
      std::default_random_engine generator;

      for (long int permutation = 0; permutation < vm["null"].as<long int>(); ++permutation)
      {
         // Generate random sample indices
         arma::uvec sampled_indices(vm["null_subsamples"].as<int>(), arma::fill::zeros);
         for (unsigned int i = 0; i < sampled_indices.n_elem; ++i)
         {
            sampled_indices(i) = i;
         }
         for (unsigned int i = sampled_indices.n_elem; i < alignment.n_rows; ++i)
         {
            std::uniform_int_distribution<int> distribution(0, i);
            unsigned int idx = distribution(generator);
            if (idx < sampled_indices.n_elem)
            {
               sampled_indices(idx) = i;
            }
         }

         sampled_indices = sort(sampled_indices);
         arma::uvec unsampled_indices(alignment.n_rows - sampled_indices.n_elem, arma::fill::zeros);
         unsigned int j = 0, k = 0;
         for (unsigned int i = 0; i < alignment.n_rows; ++i)
         {
            if (sampled_indices(j) != i)
            {
               unsampled_indices(k) = i;
               k++;
            }
            else if (j < sampled_indices.n_elem - 1) // Don't go over the end of sampled vector
            {
               j++;
            }
         }

         double D_diff = D_subsample(alignment, pairwise_distances, sampled_indices, positions, vm.count("verbose"))
                         - D_subsample(alignment, pairwise_distances, unsampled_indices, positions, vm.count("verbose"));
         std::cout << std::fixed << std::setprecision(5) << D_diff << std::endl;
      }
   }
   // Normal mode, just report D
   else
   {
      double k_hat = accu(pairwise_distances) / (num_samples * (num_samples - 1));
      double D = calc_D(k_hat, num_samples, positions, vm.count("verbose"));
      std::cout << std::fixed << std::setprecision(5) << D << std::endl;
   }

   if (vm.count("verbose"))
   {
      std::cerr << "Done.\n";
   }
}

std::tuple<long int,std::vector<int>> readCsvLine(std::string& line)
{
   std::vector<int> variant;

   std::stringstream line_stream(line);
   std::string value;

   std::getline(line_stream, value, ',');
   long int position = std::stoi(value);
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
   return std::make_tuple(position, variant);
}

double D_subsample(const arma::mat& alignment, const arma::mat& distances, const arma::uvec& sampled_indices, const std::vector<long int>& positions, const int verbose)
{
   // Take a slice from the alignment, remove any non-segregating sites
   arma::mat sub_alignment = alignment.rows(sampled_indices);
   size_t num_samples = sub_alignment.n_rows;

   arma::rowvec af = sum(sub_alignment, 0);
   int idx = 0;
   std::vector<long int> sub_pos;
   arma::uvec seg_sites(af.n_elem);
   for (unsigned int site = 0; site < af.n_elem; ++site)
   {
      if (af(site) > 0 && af(site) < sampled_indices.n_elem)
      {
         seg_sites[idx] = site;
         sub_pos.push_back(positions[site]);
         idx++;
      }
   }

   double k_hat = accu(distances(sampled_indices, sampled_indices)) / (num_samples * (num_samples - 1));

   double D = calc_D(k_hat, num_samples, sub_pos, verbose);
   return D;
}

// Pairwise distance matrix
arma::mat dist_mat(const arma::mat& alignment)
{
   arma::mat distances = arma::zeros(alignment.n_rows, alignment.n_rows);

   size_t num_samples = alignment.n_rows;
   for (size_t i = 0; i < num_samples; i++)
   {
      for (size_t j = i+1; j < num_samples; j++)
      {
         distances(i,j) = dot(alignment.row(i) - alignment.row(j), alignment.row(i) - alignment.row(j));
         distances(j,i) = distances(i,j);
      }
   }
   return distances;
}

// Calculate D from k_hat (using dist mat above) and sites
// https://en.wikipedia.org/wiki/Tajima's_D#Mathematical_details
double calc_D(const double k_hat, size_t num_samples, std::vector<long int>& positions, const int verbose)
{
   // Calculate S
   auto unique_sites = std::unique(positions.begin(), positions.end());
   double S = std::distance(positions.begin(), unique_sites);

   if (verbose)
   {
      std::cerr << "k_hat: " << k_hat << std::endl;
      std::cerr << "S: " << S << std::endl;
   }

   // Calculate D
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
   return D;
}
