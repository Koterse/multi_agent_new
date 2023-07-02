#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <numbers>
#include <random>

#include "model.h"

double a[N];
double c[N];
double Q[N][A][A];
double X[N][A];
double preX[N][A];
double action_space[M];
double action[N];
long long visit[A][A];

// multi-agent Q-learning New
int main() {
    std::mt19937 gen(19260817);
    std::binomial_distribution<int> bin(2 * S, 0.5);
    std::extreme_value_distribution<double> gumbel(
        -ETA * std::numbers::egamma_v<double>, ETA);

    // initialization stage
    initialize_Q(Q, 0);
    initialize_strategies(X);
    for (int i = 0; i < N; i++) {
        a[i] = 2;
        c[i] = 1;
    }

    initialize_action(action, 0);
    initialize_action_space(action_space);

    memcpy(preX, X, sizeof(X));

    long long t = 1;
    while (t < 1000000000LL) {
        double a_0 = 0.15;

        // policy updating
        int choice[N];
        int max_index[N];
        for (int i = 0; i < N; i++) {
            choice[i] = random_choice_from_population(X[i], gen);
            max_index[i] = -1;
            double max = -1e20;
            for (int b = 0; b < M; b++) {
                double g = gumbel(gen);
                double f = calculate_F_value(b, i, X, Q);
                if (f + g > max) {
                    max = f + g;
                    max_index[i] = b;
                }
            }
            // if (t % 100 == 0) {
            //     fprintf(stderr, "t = %lld, state = %d, i = %d, max_arg = %d
            //     max = %f\n", t, s_m, i, max_index[i], max);
            // }
        }
        for (int i = 0; i < N; ++i) {
            update_strategy(X[i], choice[i], max_index[i], calculate_nm(t));
        }

        visit[choice[0]][choice[1]]++;

        double reward[N];

        for (int i = 0; i < N; i++) {
            action[i] = action_space[choice[i]];
        }
        for (int i = 0; i < N; i++) {
            reward[i] = calculate_profit(i, t, action, a, c, a_0);
        }

        // prlong longf("%d %f\n",t,calculate_collusion_profit(t, action, a, c,
        // nash_profit, mono_profit,a_0)); if (t > T - 2000) {
        //     prlong longf("%d %f \n", t, action_space[choice[1]]);
        // }

        // updating Q value
        double tmp[2];
        double rate = calculate_learning_rate(visit[choice[0]][choice[1]]);
        for (int i = 0; i < N; i++) {
            tmp[i] = (1 - rate) * Q[i][choice[i]][choice[1 - i]] +
                     rate * (reward[i] + GAMMA * calculate_Popu(i, X, Q));
        }
        for (int i = 0; i < N; ++i) {
            Q[i][choice[i]][choice[1 - i]] = tmp[i];
        }

        if (t % 100000 == 0) {
            fprintf(stderr, "t = %lld\n", t);
            // for (int i = 0; i < N; i++) {
            //     for (int j = 0; j < S; j++) {
            //         for (int a = 0; a < A; a++) {
            //             printf("%d %f %f: %f\n", i, -0.15 + (0.15 - (-0.15))
            //             / S * j, action_space[a], X[i][j][a]);
            //         }
            //     }
            // }
            if (converge(preX, X)) {
                break;
            }
            memcpy(preX, X, sizeof(X));
        }
        if (t % 10000000 == 0) {
            printf("t = %lld\n", t);
            for (int i = 0; i < N; ++i) {
                for (int a = 0; a < A; ++a) {
                    for (int b = 0; b < A; ++b) {
                        printf("%d %d %d %lf\n", i, a, b, Q[i][a][b]);
                    }
                }
            }
            for (int i = 0; i < N; i++) {
                for (int a = 0; a < A; a++) {
                    printf("%d %f %f\n", i, action_space[a], X[i][a]);
                }
            }
        }
        fflush(stdout);
        ++t;
    }
    for (int i = 0; i < N; ++i) {
        for (int a = 0; a < A; ++a) {
            for (int b = 0; b < A; ++b) {
                printf("%d %d %d %lf\n", i, a, b, Q[i][a][b]);
            }
        }
    }
    for (int i = 0; i < N; i++) {
        for (int a = 0; a < A; a++) {
            printf("%d %f %f\n", i, action_space[a], X[i][a]);
        }
    }
    fflush(stdout);
    // for (int i = 0; i < M; i++) {
    //     for (int j = 0; j < M; j++) {
    //         action[0] = action_space[i];
    //         action[1] = action_space[j];
    //         printf("action i: %f ,action j: %f ,reward i: %f ,reward j: "
    //                "%f,average profit: % f\n ",
    //                action_space[i], action_space[j],
    //                calculate_profit(0, 0, action, a, c, 0),
    //                calculate_profit(1, 0, action, a, c, 0),
    //                (calculate_profit(0, 0, action, a, c, 0) +
    //                 calculate_profit(1, 0, action, a, c, 0)) /
    //                    2);
    //     }
    // }

    return 0;
}