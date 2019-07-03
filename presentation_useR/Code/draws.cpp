double x_prop = x_now + gsl_ran_gaussian_ziggurat(r, sigmap_beta(p));
double u = gsl_rng_uniform(r);
double x_prop = x_now + R::rnorm(0, sigmap_beta(p));
double u = R::runif(0,1);
