struct Param{
  const Rcpp::NumericVector &D;
  const Rcpp::NumericMatrix &ref;
  const Rcpp::IntegerVector &kfixed;
};

double lnl(const gsl_vector *v, void *params);

void dlnl(const gsl_vector *v, void *params, gsl_vector *df);

void lnl_dlnl(const gsl_vector *v, void *params, double *f, gsl_vector *df);

double lnl2(const gsl_vector *v, void *params);

void dlnl2(const gsl_vector *v, void *params, gsl_vector *df);

void lnl_dlnl2(const gsl_vector *v, void *params, double *f, gsl_vector *df);
