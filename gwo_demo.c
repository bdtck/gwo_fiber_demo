// gwo_fiber_demo.c
// Demo Grey Wolf Optimizer (GWO) trên bài toán khớp đường cong suy hao sợi quang
// P(L) = a * exp(-b * L) + c
// Tác giả: <Tên bạn>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N_WOLVES   20
#define DIM        3
#define MAX_ITER   200
#define N_POINTS   30

// Giới hạn cho (a, b, c)
const double LB[DIM] = {0.1,   0.001, 0.0};
const double UB[DIM] = {2.0,   0.2,   0.1};

typedef struct {
    double pos[DIM];
    double fit;
} Wolf;

// Dữ liệu sợi quang
double L_arr[N_POINTS];
double P_meas[N_POINTS];

// Tham số "thực" dùng sinh dữ liệu (để in ra so sánh)
double P0_true_g      = 1.0;
double alpha_db_true_g = 0.2;
double beta_true_g    = 0.0;
double noise_true_g   = 0.01;

// Hàm random tiện dụng
double rand01(void) {
    return (double)rand() / (double)RAND_MAX;
}

// Random Gaussian (Box-Muller) cho nhiễu, nếu bạn muốn mượt hơn
double randn(double mean, double stddev) {
    const double PI = 3.14159265358979323846;
    double u1 = rand01();
    double u2 = rand01();
    double z0 = sqrt(-2.0 * log(u1 + 1e-12)) * cos(2.0 * PI * u2);
    return mean + stddev * z0;
}

// Sinh dữ liệu mô phỏng cho sợi quang
void generate_fiber_data(void) {
    beta_true_g = alpha_db_true_g * log(10.0) / 10.0;

    for (int i = 0; i < N_POINTS; ++i) {
        L_arr[i] = 50.0 * i / (double)(N_POINTS - 1);  // 0 -> 50 km
        double P_clean = P0_true_g * exp(-beta_true_g * L_arr[i]) + noise_true_g;
        double noise = randn(0.0, 0.02);               // nhiễu ~ N(0, 0.02)
        P_meas[i] = P_clean + noise;
    }
}

// Mô hình P(L) = a * exp(-b * L) + c
double model_exp(const double params[DIM], double x) {
    double a = params[0];
    double b = params[1];
    double c = params[2];
    return a * exp(-b * x) + c;
}

// Hàm mục tiêu: MSE + ràng buộc mềm a>0, b>0, c>=0
double mse_objective(const double params[DIM]) {
    double a = params[0];
    double b = params[1];
    double c = params[2];

    // Ràng buộc đơn giản
    if (a <= 0.0 || b <= 0.0 || c < 0.0) {
        return 1e9;  // phạt lớn
    }

    double mse = 0.0;
    for (int i = 0; i < N_POINTS; ++i) {
        double y_hat = model_exp(params, L_arr[i]);
        double diff = y_hat - P_meas[i];
        mse += diff * diff;
    }
    mse /= (double)N_POINTS;
    return mse;
}

// Tìm 3 con sói tốt nhất (alpha, beta, delta)
void find_best_3(Wolf wolves[N_WOLVES], int *idx_a, int *idx_b, int *idx_d) {
    double best1 = 1e300, best2 = 1e300, best3 = 1e300;
    int i1 = -1, i2 = -1, i3 = -1;

    for (int i = 0; i < N_WOLVES; ++i) {
        double f = wolves[i].fit;
        if (f < best1) {
            best3 = best2; i3 = i2;
            best2 = best1; i2 = i1;
            best1 = f;     i1 = i;
        } else if (f < best2) {
            best3 = best2; i3 = i2;
            best2 = f;     i2 = i;
        } else if (f < best3) {
            best3 = f;     i3 = i;
        }
    }

    *idx_a = i1;
    *idx_b = i2;
    *idx_d = i3;
}

int main(void) {
    srand((unsigned int)time(NULL));

    // 1) Sinh dữ liệu
    generate_fiber_data();

    printf("=== Tham so 'thuc' dung sinh du lieu ===\n");
    printf("P0_true      = %.4f mW\n", P0_true_g);
    printf("alpha_db_true= %.4f dB/km\n", alpha_db_true_g);
    printf("beta_true    = %.6f 1/km\n", beta_true_g);
    printf("noise_floor  = %.4f mW\n\n", noise_true_g);

    // 2) Khoi tao dan soi
    Wolf wolves[N_WOLVES];
    for (int i = 0; i < N_WOLVES; ++i) {
        for (int j = 0; j < DIM; ++j) {
            double r = rand01();
            wolves[i].pos[j] = LB[j] + r * (UB[j] - LB[j]);
        }
        wolves[i].fit = mse_objective(wolves[i].pos);
    }

    // Xac dinh alpha, beta, delta
    int idx_alpha, idx_beta, idx_delta;
    find_best_3(wolves, &idx_alpha, &idx_beta, &idx_delta);

    double best_f = wolves[idx_alpha].fit;
    double best_pos[DIM];
    for (int j = 0; j < DIM; ++j) best_pos[j] = wolves[idx_alpha].pos[j];

    // 3) Vong lap GWO
    for (int iter = 1; iter <= MAX_ITER; ++iter) {
        double a = 2.0 - 2.0 * (double)iter / (double)MAX_ITER;

        // Luu lai vi tri alpha, beta, delta
        double alpha_pos[DIM], beta_pos[DIM], delta_pos[DIM];
        for (int j = 0; j < DIM; ++j) {
            alpha_pos[j] = wolves[idx_alpha].pos[j];
            beta_pos[j]  = wolves[idx_beta].pos[j];
            delta_pos[j] = wolves[idx_delta].pos[j];
        }

        // Cap nhat tung con soi
        for (int i = 0; i < N_WOLVES; ++i) {
            for (int j = 0; j < DIM; ++j) {
                double r1, r2, A1, C1, D_alpha, X1;
                double A2, C2, D_beta,  X2;
                double A3, C3, D_delta, X3;

                // Theo alpha
                r1 = rand01(); r2 = rand01();
                A1 = 2.0 * a * r1 - a;
                C1 = 2.0 * r2;
                D_alpha = fabs(C1 * alpha_pos[j] - wolves[i].pos[j]);
                X1 = alpha_pos[j] - A1 * D_alpha;

                // Theo beta
                r1 = rand01(); r2 = rand01();
                A2 = 2.0 * a * r1 - a;
                C2 = 2.0 * r2;
                D_beta = fabs(C2 * beta_pos[j] - wolves[i].pos[j]);
                X2 = beta_pos[j] - A2 * D_beta;

                // Theo delta
                r1 = rand01(); r2 = rand01();
                A3 = 2.0 * a * r1 - a;
                C3 = 2.0 * r2;
                D_delta = fabs(C3 * delta_pos[j] - wolves[i].pos[j]);
                X3 = delta_pos[j] - A3 * D_delta;

                // Vi tri moi
                double new_x = (X1 + X2 + X3) / 3.0;

                // Gioi han bien
                if (new_x < LB[j]) new_x = LB[j];
                if (new_x > UB[j]) new_x = UB[j];

                wolves[i].pos[j] = new_x;
            }

            // Tinh lai fitness
            wolves[i].fit = mse_objective(wolves[i].pos);
        }

        // Cap nhat alpha, beta, delta
        find_best_3(wolves, &idx_alpha, &idx_beta, &idx_delta);

        if (wolves[idx_alpha].fit < best_f) {
            best_f = wolves[idx_alpha].fit;
            for (int j = 0; j < DIM; ++j) best_pos[j] = wolves[idx_alpha].pos[j];
        }
    }

    // 4) In ket qua
    double a_hat = best_pos[0];
    double b_hat = best_pos[1];
    double c_hat = best_pos[2];

    double alpha_db_est = b_hat * 10.0 / log(10.0);

    printf("=== Ket qua uoc luong bang GWO (C) ===\n");
    printf("Best params (a, b, c): [%.6f, %.6f, %.6f]\n",
           a_hat, b_hat, c_hat);
    printf("  -> a_hat (P0)      ≈ %.4f mW\n", a_hat);
    printf("  -> b_hat (beta)    ≈ %.6f 1/km\n", b_hat);
    printf("  -> c_hat (floor)   ≈ %.4f mW\n", c_hat);
    printf("  -> alpha_db_est    ≈ %.4f dB/km\n", alpha_db_est);
    printf("Best MSE: %.6e\n", best_f);

    return 0;
}
