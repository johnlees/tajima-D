/*
 * File: pair.cpp
 *
 * Helper functions for the pair class
 *
 */

#include "pair.hpp"

const std::string pair_comment_default = "NA";

Pair::Pair(int number_samples)
   :_number_samples(number_samples), _bact_line(0), _human_line(0), _covars_set(0), _maf_x(0), _maf_y(0), _chisq_p(1), _lrt_p(1), _log_likelihood(0), _null_ll(0), _beta(0), _se(0), _comment(pair_comment_default), _firth(0)
{
   _x.zeros(number_samples);
   _y.zeros(number_samples);
}

// Print fields tab sep, identical to input. Doesn't print newline
std::ostream& operator<<(std::ostream &os, const Pair& p)
{
   auto maf = p.maf();
   os << std::fixed << std::setprecision(3) << p.human_line() << "\t" << p.bact_line()
      << "\t" << std::get<0>(maf) << "\t" << std::get<1>(maf)
      << "\t" << std::scientific << p.chisq_p() << "\t" << p.p_val()
      << "\t" << p.beta() << "\t" << p.comments();

   return os;
}

// Set the x and maf
void Pair::add_x(const std::vector<std::string>& variant, const long int human_line)
{
   if (variant.size() != _number_samples)
   {
      throw std::runtime_error("bacterial snps: sample size incorrect\n");
   }
   _x.zeros(variant.size());

   int i = 0;
   int missing = 0;
   for (auto it = variant.begin(); it != variant.end(); ++it)
   {
      if (*it == "0/1")
      {
         _x[i] = 1;
      }
      else if (*it == "1/1")
      {
         _x[i] = 2;
      }
      // missing as ref
      else if (*it == "./.")
      {
         missing++;
      }
      else if (*it != "0/0")
      {
         std::cerr << "none standard human snp\n";
      }
      i++;
   }

   _maf_x = (double)accu(_x)/_x.n_elem;
   _missing_x = (double)missing/_x.n_elem;

   // stats also get reset
   _human_line = human_line;
   this->reset_stats();
}

// For null ll
void Pair::add_x(const arma::mat x)
{
   _x = x;
   _human_line = 0;
}

// Set the y and maf
void Pair::add_y(const std::vector<std::string>& variant, const long int bact_line)
{
   if (variant.size() != _number_samples)
   {
      throw std::runtime_error("bacterial snps: sample size incorrect\n");
   }

   _y.zeros(variant.size());

   int i = 0;
   int missing = 0;
   for (auto it = variant.begin(); it != variant.end(); ++it)
   {
      if (*it == "1")
      {
         _y[i] = 1;
      }
      // missing as ref
      else if (*it == ".")
      {
         missing++;
      }
      else if (*it != "0")
      {
         throw std::runtime_error("none standard bacterial snp\n");
      }
      i++;
   }

   _maf_y = (double)accu(_y)/_y.n_elem;
   _missing_y = (double)missing/_y.n_elem;

   // stats also get reset
   _bact_line = bact_line;
   this->reset_stats();
}

// For null ll
void Pair::add_y(const arma::vec y)
{
   _y = y;
   _bact_line = 0;
}

// Add covariates
void Pair::add_covar(const arma::mat& covars)
{
   if (covars.n_rows != _number_samples)
   {
      throw std::runtime_error("covariates: sample size incorrect\n");
   }

   _covars = covars;
   _covars_set = 1;
}


// Add a new comment in
void Pair::add_comment(const std::string& new_comment)
{
   if (_comment == pair_comment_default)
   {
      _comment = new_comment;
   }
   else
   {
      _comment += "," + new_comment;
   }
}

// Get covars
arma::mat Pair::get_covars()
{
   if (!(_covars_set))
   {
      throw std::logic_error("Tried to access pair covars when they have not been set");
   }

   return _covars;
}

// Get the design matrix
arma::mat Pair::get_x_design()
{
   arma::mat x_design = arma::join_rows(arma::ones<arma::vec>(_number_samples), _x);
   if (_covars_set)
   {
      x_design = arma::join_rows(x_design, _covars);
   }

   return x_design;
}

void Pair::reset_stats()
{
   // null_ll is retained
   _chisq_p = 1;
   _lrt_p = 1;
   _log_likelihood = 0;
   _beta = 0;
   _se = 0;
   _comment = pair_comment_default;
   _firth = 0;
}

