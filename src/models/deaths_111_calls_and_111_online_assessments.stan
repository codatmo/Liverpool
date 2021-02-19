functions {
  real[] seeiittd(real time,
                  real[] state,
                  real[] params,
                  real[] real_data,
                  int[] integer_data) {

    // Unpack integer data values
    int T = integer_data[1];
    int n_beta_pieces = integer_data[2];
    int n_disease_states = integer_data[3];

    // Unpack real data values
    real beta_left_t[n_beta_pieces] = real_data[1:n_beta_pieces];
    real beta_right_t[n_beta_pieces] = real_data[n_beta_pieces+1:2*n_beta_pieces];
    real population = real_data[2*n_beta_pieces+1];

    // Unpack parameter values
    real beta_left[n_beta_pieces] = params[1:n_beta_pieces];
    real grad_beta[n_beta_pieces] = params[n_beta_pieces+1:2*n_beta_pieces];
    real nu = params[2*n_beta_pieces+1];
    real gamma = params[2*n_beta_pieces+2];
    real kappa = params[2*n_beta_pieces+3];
    real omega = params[2*n_beta_pieces+4];

    // Unpack state
    real S = state[1];
    real E1 = state[2];
    real E2 = state[3];
    real I1 = state[4];
    real I2 = state[5];
    real T1 = state[6];
    real T2 = state[7];
    real D = state[8];

    real infection_rate;

    for (i in 1:n_beta_pieces) {
      if(time >= beta_left_t[i] && time < beta_right_t[i]) {
        real beta = grad_beta[i] * (time - beta_left_t[i]) + beta_left[i];
        infection_rate = beta * (I1 + I2) * S / population;
      }
    }

    real nuE1 = nu * E1;
    real nuE2 = nu * E2;
    real gammaI1 = gamma * I1;
    real gammaI2 = gamma * I2;
    real kappaT1 = kappa * T1;
    real kappaT2 = kappa * T2;

    real dS_dt = -infection_rate;
    real dE1_dt = infection_rate - nuE1;
    real dE2_dt = nuE1 - nuE2;
    real dI1_dt = nuE2 - gammaI1;
    real dI2_dt = gammaI1 - gammaI2;
    real dT1_dt = gammaI2 * omega - kappaT1;
    real dT2_dt = kappaT1 - kappaT2;
    real dD_dt = kappaT2;

    return {dS_dt, dE1_dt, dE2_dt, dI1_dt, dI2_dt, dT1_dt, dT2_dt, dD_dt};
  }

  real[ , ] integrate_ode_explicit_trapezoidal(real[] initial_state, real initial_time, real[] times, real[] params, real[] real_data, int[] integer_data) {
    real h;
    vector[size(initial_state)] dstate_dt_initial_time;
    vector[size(initial_state)] dstate_dt_tidx;
    vector[size(initial_state)] k;
    real state_estimate[size(times),size(initial_state)];

    h = times[1] - initial_time;
    dstate_dt_initial_time = to_vector(seeiittd(initial_time, initial_state, params, real_data, integer_data));
    k = h*dstate_dt_initial_time;
    state_estimate[1,] = to_array_1d(to_vector(initial_state) + h*(dstate_dt_initial_time + to_vector(seeiittd(times[1], to_array_1d(to_vector(initial_state)+k), params, real_data, integer_data)))/2);

    for (tidx in 1:size(times)-1) {
      h = (times[tidx+1] - times[tidx]);
      dstate_dt_tidx = to_vector(seeiittd(times[tidx], state_estimate[tidx], params, real_data, integer_data));
      k = h*dstate_dt_tidx;
      state_estimate[tidx+1,] = to_array_1d(to_vector(state_estimate[tidx,]) + h*(dstate_dt_tidx + to_vector(seeiittd(times[tidx+1], to_array_1d(to_vector(state_estimate[tidx,])+k), params, real_data, integer_data)))/2);
    }

    return state_estimate;
  }
}
data {
  real initial_time;
  int<lower=1> n_beta_pieces;
  real<lower=0> beta_left_t[n_beta_pieces];
  real<lower=0> beta_right_t[n_beta_pieces];
  int<lower=1> n_rho_calls_111_pieces;
  int<lower=0> rho_calls_111_left_t[n_rho_calls_111_pieces];
  int<lower=0> rho_calls_111_right_t[n_rho_calls_111_pieces];
  int<lower=1> n_rho_online_assessments_111_pieces;
  int<lower=0> rho_online_assessments_111_left_t[n_rho_online_assessments_111_pieces];
  int<lower=0> rho_online_assessments_111_right_t[n_rho_online_assessments_111_pieces];
  int<lower=1> T;
  real times[T];
  int<lower=1> n_disease_states;
  real<lower=0> population;
  int<lower=1> deaths_length;
  int<lower=1> deaths_starts[deaths_length];
  int<lower=1> deaths_stops[deaths_length];
  int<lower=0> deaths[deaths_length];
  int<lower=1> calls_111_length;
  int<lower=1> calls_111_start;
  int<lower=0> calls_111[calls_111_length];
  int<lower=1> online_assessments_111_length;
  int<lower=1> online_assessments_111_start;
  int<lower=0> online_assessments_111[online_assessments_111_length];
  int real_data_length;
  real real_data[real_data_length];
  int integer_data_length;
  int integer_data[integer_data_length];
}
transformed data {
  real mu_dL = 4.00;
  real sigma_dL = 0.20;
  real mu_dI = 3.06;
  real sigma_dI = 0.21;
  real mu_dT = 16.00;
  real sigma_dT = 0.71;
  int max_lag = 13;
}
parameters {
  real<lower=0, upper=1> initial_state_raw[2];
  real<lower=0> beta_left[n_beta_pieces];
  real<lower=0> beta_right[n_beta_pieces];
  real<lower=0> dL;
  real<lower=0> dI;
  real<lower=0> dT;
  real<lower=0, upper=1> omega;
  real<lower=0> reciprocal_phi_deaths;
  real<lower=0> reciprocal_phi_calls_111;
  real<lower=0> reciprocal_phi_online_assessments_111;
  real<lower=0> rho_calls_111[n_rho_calls_111_pieces];
  real<lower=0> rho_online_assessments_111[n_rho_online_assessments_111_pieces];
  simplex[max_lag+1] lag_weights_calls_111;
  simplex[max_lag+1] lag_weights_online_assessments_111;
}
transformed parameters {
  real initial_state[n_disease_states];
  real grad_beta[n_beta_pieces];
  real nu;
  real gamma;
  real kappa;
  real phi_deaths;
  real phi_calls_111;
  real phi_online_assessments_111;
  real state_estimate[T,n_disease_states];
  vector[T+1] S;
  vector[T+1] E1;
  vector[T+1] E2;
  vector[T+1] I1;
  vector[T+1] I2;
  vector[T+1] T1;
  vector[T+1] T2;
  vector[T+1] D;
  vector[T] daily_infections;
  vector[T] daily_deaths;
  vector[T] calls_111_lagged_daily_infections;
  vector[T] online_assessments_111_lagged_daily_infections;
  vector[T] daily_calls_111;
  vector[T] daily_online_assessments_111;

  initial_state[1] = (population-5.0)*initial_state_raw[1] + 1.0;
  initial_state[2] = (population-5.0)*(1.0-initial_state_raw[1])*initial_state_raw[2]/2.0 + 1.0;
  initial_state[3] = (population-5.0)*(1.0-initial_state_raw[1])*initial_state_raw[2]/2.0 + 1.0;
  initial_state[4] = (population-5.0)*(1.0-initial_state_raw[1])*(1.0-initial_state_raw[2])/2.0 + 1.0;
  initial_state[5] = (population-5.0)*(1.0-initial_state_raw[1])*(1.0-initial_state_raw[2])/2.0 + 1.0;
  initial_state[6] = 0.0;
  initial_state[7] = 0.0;
  initial_state[8] = 0.0;
  grad_beta = to_array_1d((to_vector(beta_right) - to_vector(beta_left))./(to_vector(beta_right_t) - to_vector(beta_left_t)));
  nu = 2.0/dL;
  gamma = 2.0/dI;
  kappa = 2.0/dT;
  phi_deaths = 1.0 / reciprocal_phi_deaths;
  phi_calls_111 = 1.0 / reciprocal_phi_calls_111;
  phi_online_assessments_111 = 1.0 / reciprocal_phi_online_assessments_111;

  {
    real params[2*n_beta_pieces+4];
    params[1:n_beta_pieces] = beta_left;
    params[n_beta_pieces+1:2*n_beta_pieces] = grad_beta;
    params[2*n_beta_pieces+1] = nu;
    params[2*n_beta_pieces+2] = gamma;
    params[2*n_beta_pieces+3] = kappa;
    params[2*n_beta_pieces+4] = omega;

    state_estimate = integrate_ode_explicit_trapezoidal(initial_state, initial_time, times, params, real_data, integer_data);
  }

  S = append_row(initial_state[1], to_vector(state_estimate[, 1]));
  E1 = append_row(initial_state[2], to_vector(state_estimate[, 2]));
  E2 = append_row(initial_state[3], to_vector(state_estimate[, 3]));
  I1 = append_row(initial_state[4], to_vector(state_estimate[, 4]));
  I2 = append_row(initial_state[5], to_vector(state_estimate[, 5]));
  T1 = append_row(initial_state[6], to_vector(state_estimate[, 6]));
  T2 = append_row(initial_state[7], to_vector(state_estimate[, 7]));
  D = append_row(initial_state[8], to_vector(state_estimate[, 8]));

  daily_infections = S[:T] - S[2:] + machine_precision();
  daily_deaths = D[2:] - D[:T];

  calls_111_lagged_daily_infections = lag_weights_calls_111[1]*daily_infections;

  for (i in 1:max_lag) {
    calls_111_lagged_daily_infections += lag_weights_calls_111[i+1]*append_row(rep_vector(0.0, i), daily_infections[:T-i]);
  }

  daily_calls_111 = rep_vector(0.0, T);

  for (i in 1:n_rho_calls_111_pieces) {
    daily_calls_111[rho_calls_111_left_t[i]:rho_calls_111_right_t[i]-1] = calls_111_lagged_daily_infections[rho_calls_111_left_t[i]:rho_calls_111_right_t[i]-1] * rho_calls_111[i];
  }

  online_assessments_111_lagged_daily_infections = lag_weights_online_assessments_111[1]*daily_infections;

  for (i in 1:max_lag) {
    online_assessments_111_lagged_daily_infections += lag_weights_online_assessments_111[i+1]*append_row(rep_vector(0.0, i), daily_infections[:T-i]);
  }

  daily_online_assessments_111 = rep_vector(0.0, T);

  for (i in 1:n_rho_online_assessments_111_pieces) {
    daily_online_assessments_111[rho_online_assessments_111_left_t[i]:rho_online_assessments_111_right_t[i]-1] = online_assessments_111_lagged_daily_infections[rho_online_assessments_111_left_t[i]:rho_online_assessments_111_right_t[i]-1] * rho_online_assessments_111[i];
  }
}
model {
  initial_state_raw[1] ~ beta(5.0, 0.5);
  initial_state_raw[2] ~ beta(1.1, 1.1);
  beta_left ~ normal(0, 0.5);
  beta_right ~ normal(0, 0.5);
  dL ~ normal(mu_dL, sigma_dL);
  dI ~ normal(mu_dI, sigma_dI);
  dT ~ normal(mu_dT, sigma_dT);
  omega ~ beta(100, 9803);
  reciprocal_phi_deaths ~ exponential(5);
  reciprocal_phi_calls_111 ~ exponential(5);
  reciprocal_phi_online_assessments_111 ~ exponential(5);
  rho_calls_111 ~ normal(0, 0.5);
  rho_online_assessments_111 ~ normal(0, 0.5);
  lag_weights_calls_111 ~ dirichlet(rep_vector(0.1, max_lag+1));
  lag_weights_online_assessments_111 ~ dirichlet(rep_vector(0.1, max_lag+1));

  for (i in 1:deaths_length) {
    target += 7.0*neg_binomial_2_lpmf(deaths[i] | sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]), phi_deaths);
  }

  for (i in 1:n_rho_calls_111_pieces) {
    target += neg_binomial_2_lpmf(calls_111[rho_calls_111_left_t[i]-calls_111_start+1:rho_calls_111_right_t[i]-calls_111_start] | daily_calls_111[rho_calls_111_left_t[i]:rho_calls_111_right_t[i]-1], phi_calls_111);
  }

  for (i in 1:n_rho_online_assessments_111_pieces) {
    target += neg_binomial_2_lpmf(online_assessments_111[rho_online_assessments_111_left_t[i]-online_assessments_111_start+1:rho_online_assessments_111_right_t[i]-online_assessments_111_start] | daily_online_assessments_111[rho_online_assessments_111_left_t[i]:rho_online_assessments_111_right_t[i]-1], phi_online_assessments_111);
  }
}
generated quantities {
  vector[T-1] growth_rate = (log(daily_infections[2:]) - log(daily_infections[:T-1]))*100;

  int pred_deaths[deaths_length];

  for (i in 1:deaths_length) {
    pred_deaths[i] = neg_binomial_2_rng(sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]), phi_deaths);
  }

  int pred_calls_111[calls_111_length];

  for (i in 1:n_rho_calls_111_pieces) {
    pred_calls_111[rho_calls_111_left_t[i]-calls_111_start+1:rho_calls_111_right_t[i]-calls_111_start] = neg_binomial_2_rng(daily_calls_111[rho_calls_111_left_t[i]:rho_calls_111_right_t[i]-1], phi_calls_111);
  }

  int pred_online_assessments_111[online_assessments_111_length];

  for (i in 1:n_rho_online_assessments_111_pieces) {
    pred_online_assessments_111[rho_online_assessments_111_left_t[i]-online_assessments_111_start+1:rho_online_assessments_111_right_t[i]-online_assessments_111_start] = neg_binomial_2_rng(daily_online_assessments_111[rho_online_assessments_111_left_t[i]:rho_online_assessments_111_right_t[i]-1], phi_online_assessments_111);
  }
}
