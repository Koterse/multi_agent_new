#include <iostream>
#include <cstring>
#include "model.h"

double reward[A][A][A][3];
double action_space[A];
double p[3][A];

void expected(double p[3][A], double res[3]) {
    res[0] = res[1] = res[2] = 0;
    for (int i = 0; i < A; ++i) {
        for (int j = 0; j < A; ++j) {
            for (int k = 0; k < A; ++k) {
                double pp = p[0][i] * p[1][j] * p[2][k];
                res[0] += pp * reward[i][j][k][0];
                res[1] += pp * reward[i][j][k][1];
                res[2] += pp * reward[i][j][k][2];
            }
        }
    }
}

int main() {
    initialize_action_space(action_space);
    double a[3] = {2, 2, 2};
    double c[3] = {1, 1, 1};
    double a0 = 0;
    for (int i = 0; i < A; ++i) {
        for (int j = 0; j < A; ++j) {
            for (int k = 0; k < A; ++k) {
                double actions[3] = {action_space[i], action_space[j], action_space[k]};
                reward[i][j][k][0] = calculate_profit(0, 0, actions, a, c, a0);
                reward[i][j][k][1] = calculate_profit(1, 0, actions, a, c, a0);
                reward[i][j][k][2] = calculate_profit(2, 0, actions, a, c, a0);
            }
        }
    }
    // printf("%f\n", reward[8][8][8][0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < A; ++j) {
            int x;
            double y, z;
            scanf("%d %lf %lf", &x, &y, &z);
            p[i][j] = z;
            // printf("***%d %d %f %f\n", i, j, z, p[i][j]);
        }
    }
    // printf("%f\n", p[0][2]);
    double now[3];
    expected(p, now);
    for (int i = 0; i < 3; ++i) {
        printf("%f ", now[i]);
    }
    puts("");
    for (int i = 0; i < 3; ++i) {
        double q[3][A];
        double eps = 0;
        memcpy(q, p, sizeof(q));
        for (int j = 0; j < A; ++j) {
            memset(q[i], 0, sizeof(q[i]));
            q[i][j] = 1;
            double then[3];
            expected(q, then);
            if (then[i] > now[i]) {
                eps = std::max(eps, (then[i] - now[i]) / now[i]);
            }
        }
        printf("%d %lf\n", i, eps);
    }
    return 0;
}
