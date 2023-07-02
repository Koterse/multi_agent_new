#ifndef MODEL_H
#define MODEL_H

#include <random>

// #define A_0 0
#define MU 0.25
#define N 2
#define DELTA 0.95
#define XI 0.2
#define M 15 //size of action space
#define N_0 1000
#define T 1e8
#define P_N 1.473
#define P_M (1 + 0.473 * 2)
#define A 15
#define S 6
#define BETA 0.1
#define ALPHA 0.1
#define GAMMA 0.1
#define KAPPA 1e4
#define ETA 0.01
#define EPSILON 5*1e-3L

long long calculate_power(long long base, long long exponent);

double calculate_demand(int item, long long time, double action[N], double *a, double a_0);

double calculate_profit(int item, long long time, double action[N], double *a, double *c, double a_0);

double calculate_Bertrand_Nash();

double calculate_monopoly();

void initialize_Q(double Q[N][A][A], double initial_value);

void initialize_strategies(double X[N][A]);

void initialize_action(double action[N], double initial_value);

void initialize_action_space(double *action_space);

void find_state_price(int index, double *result, double *action_space);

int random_choice_from_population(double X[A], std::mt19937 &gen);

double calculate_F_value(int b, int item, double X[N][A], double Q[N][A][A]);

void update_strategy(double X[A], int choice, int optimal_choice, double nm);

double calculate_Popu(int item, double X[N][A], double Q[N][A][A]);

double calculate_collusion_profit(long long time, double action[N], double *a, double *c, double nash_profit, double monopoly_profit, double a_0);

double calculate_coprofit(double *a, double *c, double p_n, double a_0);

double calculate_nm(long long n);
double calculate_learning_rate(long long vis);

bool converge(double X[N][A], double Y[N][A]);

#endif