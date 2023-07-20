#include <iostream>

#include "model.h"

double reward[A][A][A][3];
double action_space[A];

bool is_nash(int i, int j, int k) {
    for (int a = 0; a < A; ++a) {
        if (reward[a][j][k][0] > reward[i][j][k][0]) {
            return false;
        }
        if (reward[i][a][k][1] > reward[i][j][k][1]) {
            return false;
        }
        if (reward[i][j][a][2] > reward[i][j][k][2]) {
            return false;
        }
    }
    return true;
}

int main() {
    initialize_action_space(action_space);
    // for (int i = 0; i < A; ++i) {
    //     printf("%lf\n", action_space[i]);
    // }
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
    for (int i = 0; i < A; ++i) {
        for (int j = 0; j < A; ++j) {
            for (int k = 0; k < A; ++k) {
                if (is_nash(i, j, k)) {
                    printf(
                        "%lf %lf %lf -> %lf %lf %lf \n",
                        action_space[i],
                        action_space[j],
                        action_space[k],
                        reward[i][j][k][0],
                        reward[i][j][k][1],
                        reward[i][j][k][2]
                    );
                }
            }
        }
    }
    return 0;
}
