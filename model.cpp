#include "model.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <random>

long long calculate_power(long long base, long long exponent) {
    long long result = 1;
    for (long long i = 0; i < exponent; i++) {
        result *= base;
    }
    return result;
}

double calculate_demand(int item, long long time, double action[N], double *a, double a_0) {
    double sum = 0;
    for (int j = 0; j < N; j++) {
        sum += exp((a[j] - action[j]) / MU);
    }
    return exp((a[item] - action[item]) / MU) / (sum + exp(a_0 / MU));
}

double
calculate_profit(int item, long long time, double action[N], double *a, double *c, double a_0) {
    return (action[item] - c[item]) * calculate_demand(item, time, action, a, a_0);
}

double calculate_Bertrand_Nash() {
    double Bertrand_Nash = 0;
    // to be implemented using foc
    return Bertrand_Nash;
}

double calculate_monopoly() {
    double monopoly = 0;
    // to be implemented using foc_monopoly
    return monopoly;
}


void initialize_Q(double Q[N][A][A][A], double initial_value) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < A; j++) {
            for (int k = 0; k < A; k++) {
                for (int o = 0; o < A; o++) {
                    Q[i][j][k][o] = initial_value;
                }
            }
        }
    }
}


void initialize_strategies(double X[N][A]) {
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < A; k++) {
            X[i][k] = 1.0 / A;
        }
    }
}

void initialize_action(double action[N], double initial_value) {
    for (int i = 0; i < N; i++) {
        action[i] = initial_value;
    }
}


// according to assumption 4.1

double calculate_nm(long long t) {
    return 1e-5 * pow(t, 0.6) + N_0;
    // return N_0;
}

double calculate_learning_rate(long long vis) {
    return 1 / (1e-5 * pow(vis + 1, 0.6) + 10);
    // return 1e-1;
}

void initialize_action_space(double *action_space) {
    double lower = P_N - (P_M - P_N) / (M - 3);
    double upper = P_M + (P_M - P_N) / (M - 3);
    // double lower = 1;
    // double upper = 2;
    for (long long i = 0; i < M; i++) {
        action_space[i] = lower + i * ((upper - lower) / (M - 1));
    }
}
// only applicable to two agent learg
// this state is defined using memory
void find_state_price(int index, double *result, double *action_space) {
    int a = index / M;
    int b = index % M;
    result[0] = action_space[a];
    result[1] = action_space[b];
}

int random_choice_from_population(double X[A], std::mt19937 &gen) {
    std::discrete_distribution<int> distrib(X, X + A);
    return distrib(gen);
}

// only applicable to two-agents learning
double calculate_F_value(int b, int i, double X[N][A], double Q[N][A][A][A]) {
    double F = 0;
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < M; ++k) {
            F += X[(i + 1) % N][j] * X[(i + 2) % N][j] * Q[i][b][j][k];
        }
    }
    return F;
}

void update_strategy(double X[A], int choice, int optimal_choice, double nm) {
    X[optimal_choice] += 1 / nm;
    X[choice] -= 1 / nm;
}

// only applicable to two-agent learning
double calculate_Popu(int item, double X[N][A], double Q[N][A][A][A]) {
    double Popu = 0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < M; ++k) {
                // if(std::isnan(X[0][next_state][i])) {
                //     fprlong longf(stderr, "wow\n");
                // }
                // if(std::isnan(X[1][next_state][i])) {
                //     fprintf(stderr, "wow\n");
                // }
                Popu += X[item][i] * X[(item + 1) % N][j] * X[(item + 2) % N][k] * Q[item][i][j][k];
                // if (std::isnan(Popu)) {
                //     fprintf(stderr, "wow %lf %lf %lf\n", X[0][next_state][i],
                //     X[1][next_state][j], Q[item][next_state][i][j]); return Popu;
                // }
            }
        }
    }
    return Popu;
}

double calculate_collusion_profit(
    long long time,
    double action[N],
    double *a,
    double *c,
    double nash_profit,
    double monopoly_profit,
    double a_0
) {
    double triangle = 0;
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += calculate_profit(i, time, action, a, c, a_0);
    }
    double average_profit = sum / N;
    triangle = (average_profit - nash_profit) / (monopoly_profit - nash_profit);
    return triangle;
}

double calculate_coprofit(double *a, double *c, double p_n, double a_0) {
    double sum = 0;
    for (int j = 0; j < N; j++) {
        sum += exp((a[j] - p_n) / MU);
    }
    // printf("%f\n", p_n - c[0]);
    // printf("%f\n", exp(a[0] - p_n) / (sum + exp(A_0 / MU)));
    return exp((a[0] - p_n) / MU) / (sum + exp(a_0 / MU)) * (p_n - c[0]);
}

bool converge(double X[N][A], double Y[N][A]) {
    for (int i = 0; i < N; ++i) {
        long double dot = 0;
        for (int k = 0; k < A; ++k) {
            dot += powl(X[i][k] - Y[i][k], 2);
        }
        if (dot >= EPSILON * EPSILON) {
            return false;
        }
    }
    return true;
}
