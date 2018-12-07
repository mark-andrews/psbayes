data {
  int<lower=0> n;
  int<lower=0> m;
}

parameters {
  real<lower=0, upper=1> theta;
}

model {
  theta ~ beta(1, 1);
  m ~ binomial(n, theta);
}

