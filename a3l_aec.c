/*
File:       a3l_aec.c
Author:     Aironi Carlo
Date:       08/04/2025

Implementation of the (AU)MDF algorithm for acoustic echo cancellation, described in:

- J. S. Soo, K. K. Pang Multidelay block frequency adaptive filter,  IEEE Trans. Acoust.
Speech Signal Process., Vol. ASSP-38, No. 2, February 1990.

Robustness to double-talk is achieved using a variable learning rate as described in:

- Valin, J.-M., On Adjusting the Learning Rate in Frequency Domain Echo 
Cancellation With Double-Talk. IEEE Transactions on Audio, Speech and Language Processing,
Vol. 15, No. 3, pp. 1030-1034, 2007. http://people.xiph.org/~jm/papers/valin_taslp2006.pdf

Floating point 32bit (float) version - Il programma usa esclusivamente operazioni a precisione singola (float)
senza conversioni implicite a double.

*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "a3l_aec.h"

#ifndef M_PI
#define M_PI 3.141592f
#endif

#define MDF                                                                   // define MDF or AUMDF

// Allocazione e inizializzazione della struttura
void a3l_aec_state_init(a3l_aec_state *st, int frame_size, int filter_length, int sample_rate) {
    // Inizializzazione degli scalari
    st->frame_size = frame_size;                // typical value: 128
    st->window_size = 2 * frame_size;           // typical value: 256
    st->M = filter_length / frame_size;         // typical value: 16
    st->update_foreground = 0;
    st->reset_background = 0;
    st->cancel_count = 0;
    st->sum_adapt = 0.0f;
    st->saturated = 0;
    st->screwed_up = 0;
    st->spec_average = (float)st->frame_size / sample_rate;
    st->beta0 = (2.0f * st->frame_size) / sample_rate;
    st->beta_max = (0.5f * st->frame_size) / sample_rate;
    st->leak_estimate = 0.0f;
    st->prop_sum = 0.0f;

    // Allocazione e inizializzazione degli array 1D (calloc alloca e inizializza a zero)
    int N = st->window_size;
    int M = st->M;
    st->e = (float*)calloc(N, sizeof(float));
    st->x = (float*)calloc(N, sizeof(float));
    st->ifft_buf = (float*)calloc(N, sizeof(float));
    st->input = (float*)calloc(st->frame_size, sizeof(float));
    st->out = (float*)calloc(st->frame_size, sizeof(float));
    st->y = (float*)calloc(N, sizeof(float));
    st->last_y = (float*)calloc(N, sizeof(float));
    st->Yf = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->Rf = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->Xf = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->Yh = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->Eh = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->Eh_cur = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->Yh_cur = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->fft_buf = (float _Complex*)calloc(N, sizeof(float _Complex));
    st->Y = (float _Complex*)calloc(N, sizeof(float _Complex));
    st->E = (float _Complex*)calloc(N, sizeof(float _Complex));
    st->PHI = (float _Complex*)calloc(N, sizeof(float _Complex));
    st->power = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->power_1 = (float*)calloc(st->frame_size + 1, sizeof(float));
    st->window = (float*)calloc(N, sizeof(float));
    st->prop = (float*)calloc(M, sizeof(float));
    st->max_sum = (float*)calloc(M, sizeof(float));
    st->wtmp = (float*)calloc(N, sizeof(float));

    for (int i = 0; i < st->frame_size + 1; i++) {                              // power_1 è inizializzato a ones
        st->power_1[i] = 1.0f;
    }

    // Allocazione e inizializzazione degli array 2D a valori complessi
    st->X = (float _Complex**)calloc(M + 1, sizeof(float _Complex*));         // in C, X è [M+1, N], in MATLAB è [N, M+1]
    for (int i = 0; i < M + 1; i++) {
        st->X[i] = (float _Complex*)calloc(N, sizeof(float _Complex));
    }

    st->W = (float _Complex**)calloc(M, sizeof(float _Complex*));             // in C, W è [M, N], in MATLAB è [N, M]
    for (int i = 0; i < M; i++) {
        st->W[i] = (float _Complex*)calloc(N, sizeof(float _Complex));
    }

    st->foreground = (float _Complex**)calloc(M, sizeof(float _Complex*));    // in C, foreground è [M, N], in MATLAB è [N, M]
    for (int i = 0; i < M; i++) {
        st->foreground[i] = (float _Complex*)calloc(N, sizeof(float _Complex));
    }

    // Inizializzazione della finestra con il profilo di Hanning
    for (int i = 0; i < N; i++) {
        st->window[i] = 0.5f - 0.5f * cosf(2.0f * M_PI * ((float)(i + 1) - 1) / N);
    }

    // Altre inizializzazioni
    float decay = exp(-2.4f / M);
    st->prop[0] = 0.7f;
    for (int i = 1; i < M; i++) {
        st->prop[i] = st->prop[i - 1] * decay;
    }
    float sum_prop = 0.0f;
    for (int i = 0; i < M; i++) {
        sum_prop += st->prop[i];
    }
    for (int i = 0; i < M; i++) {
        st->prop[i] = 0.8f * st->prop[i] / sum_prop;
    }

    // Altri parametri
    st->adapted = 0;
    st->Pey = 1.0f;
    st->Pyy = 1.0f;
    st->Davg1 = 0.0f;
    st->Davg2 = 0.0f;
    st->Dvar1 = 0.0f;
    st->Dvar2 = 0.0f;
}

// Reset della struttura
void a3l_aec_state_reset(a3l_aec_state *st) {
    int N = st->window_size;
    int M = st->M;
    // Inizializzazione degli scalari
    st->update_foreground = 0;
    st->reset_background = 0;
    st->cancel_count = 0;
    st->screwed_up = 0;
    st->saturated = 0;
    st->sum_adapt = 0.0f;
    st->leak_estimate = 0.0f;
    st->prop_sum = 0.0f;
    // Reset degli array 1D
    memset(st->e, 0.0f, N * sizeof(float));
    memset(st->x, 0.0f, N * sizeof(float));
    memset(st->input, 0.0f, st->frame_size * sizeof(float));
    memset(st->out, 0.0f, st->frame_size * sizeof(float));
    memset(st->y, 0.0f, N * sizeof(float));
    memset(st->last_y, 0.0f, N * sizeof(float));
    memset(st->Yf, 0.0f, (st->frame_size + 1) * sizeof(float));
    memset(st->Rf, 0.0f, (st->frame_size + 1) * sizeof(float));
    memset(st->Xf, 0.0f, (st->frame_size + 1) * sizeof(float));
    memset(st->Yh, 0.0f, (st->frame_size + 1) * sizeof(float));
    memset(st->Eh, 0.0f, (st->frame_size + 1) * sizeof(float));
    memset(st->Y, 0.0f, N * sizeof(float _Complex));
    memset(st->E, 0.0f, N * sizeof(float _Complex));
    memset(st->PHI, 0.0f, N * sizeof(float _Complex));
    memset(st->power, 0.0f, (st->frame_size + 1) * sizeof(float));
    memset(st->power_1, 0.0f, (st->frame_size + 1) * sizeof(float));              // un po' inutile perché dopo lo riempio con 1...
    memset(st->Eh_cur, 0.0f, (st->frame_size + 1) * sizeof(float));
    memset(st->Yh_cur, 0.0f, (st->frame_size + 1) * sizeof(float));
    memset(st->prop, 0.0f, M * sizeof(float));
    memset(st->max_sum, 0.0f, M * sizeof(float));
    memset(st->wtmp, 0.0f, N * sizeof(float));

    for (int i = 0; i < st->frame_size + 1; i++) {                              // power_1 era inizializzato a ones
        st->power_1[i] = 1.0f;
    }
    // Reset degli array 2D a valori complessi
    for (int i = 0; i < M + 1; i++) {
        memset(st->X[i], 0.0f, N * sizeof(float _Complex));
    }

    for (int i = 0; i < M; i++) {
        memset(st->W[i], 0.0f, N * sizeof(float _Complex));
    }

    for (int i = 0; i < M; i++) {
        memset(st->foreground[i], 0.0f, N * sizeof(float _Complex));
    }
    
    float decay = exp(-2.4f / M);
    st->prop[0] = 0.7f;
    for (int i = 1; i < M; i++) {
        st->prop[i] = st->prop[i - 1] * decay;
    }
    float sum_prop = 0.0f;
    for (int i = 0; i < M; i++) {
        sum_prop += st->prop[i];
    }
    for (int i = 0; i < M; i++) {
        st->prop[i] = 0.8f * st->prop[i] / sum_prop;
    }
    
    st->adapted = 0;
    st->Pey = 1.0f;
    st->Pyy = 1.0f;
    st->Davg1 = 0.0f;
    st->Davg2 = 0.0f;
    st->Dvar1 = 0.0f;
    st->Dvar2 = 0.0f;
    
}

// Deallocazione della struttura
void a3l_aec_state_free(a3l_aec_state *st) {
    // Dealloca la memoria per ciascun array membro
    free(st->e);
    free(st->x);
    free(st->input);
    free(st->out);
    free(st->y);
    free(st->last_y);
    free(st->Yf);
    free(st->Rf);
    free(st->Xf);
    free(st->Yh);
    free(st->Eh);
    free(st->Eh_cur);
    free(st->Yh_cur);
    for (int i = 0; i < st->M + 1; i++) {
        free(st->X[i]);                             // Libera la memoria per ogni riga
    }
    free(st->X);                                    // Libera la memoria per il puntatore principale
    for (int i = 0; i < st->M; i++) {
        free(st->W[i]);
        free(st->foreground[i]);
    }
    free(st->W);
    free(st->foreground);
    free(st->Y);
    free(st->E);
    free(st->fft_buf);
    free(st->ifft_buf);
    free(st->PHI);
    free(st->power);
    free(st->power_1);
    free(st->window);
    free(st->prop);
    free(st->wtmp);
    free(st->max_sum);
    // Dopo aver liberato i membri dinamici, libera la struttura stessa
    free(st);
}

// conversione da int16_t a float e scalatura da [-32768,32767] a [-1,1]
static void int2float(int16_t *input, float *output, size_t numSamples) {
    for (size_t i = 0; i < numSamples; ++i) {
        output[i] = (float)input[i] / 32768.0f;
    }
}

// conversione da float a int16_t e scalatura da [-1,1] a [-32768,32767]
static void float2int(float *input, int16_t *output, size_t numSamples) {
    for (size_t i = 0; i < numSamples; ++i) {
        output[i] = (int16_t)round(input[i] * 32768.0f);
        if (output[i] > 32767) {
            output[i] = 32767;
        } else if (output[i] < -32768) {
            output[i] = -32768;
        }
    }
}

// Funzione circshift per matrici complesse
void circshift(a3l_aec_state* st) {
    // Alloca memoria per una riga temporanea
    float _Complex* temp = (float _Complex*)calloc(st->window_size, sizeof(float _Complex));
    // Salva l'ultima riga in `temp`
    for (int i = 0; i < st->window_size; i++) {
        temp[i] = st->X[st->M][i];
    }
    // Sposta ogni colonna in basso
    for (int j = st->M; j > 0; j--) {
        for (int i = 0; i < st->window_size; i++) {
            st->X[j][i] = st->X[j - 1][i];
        }
    }
    // Copia `temp` nella prima riga
    for (int i = 0; i < st->window_size; i++) {
        st->X[0][i] = temp[i];
    }
    // Libera la memoria della colonna temporanea
    free(temp);
}

// Funzione che rende coniugato simmetrico un array complesso di lunghezza L (pari)
static void conj_sym(float _Complex *Arr, int L) {
    // L è la lunghezza totale dell'array (uguale alla lunghezza del segnale originale)
    // Le prime 129 posizioni contengono i campioni dello spettro (frequenze positive)
    // Copia e coniuga la parte del secondo spettro (frequenze negative)
    for (int i = 1; i <= L/2; i++) {
        Arr[L-i] = conjf(Arr[i]);
    }
}

// Funzione che moltiplica per una costante tutti gli elementi di un array a valori complessi di lunghezza L
// K1 = numeratore
// K2 = denominatore
static void arr_scaling_cplx(float _Complex *Arr, int L, float K1, float K2) {
    for (int i = 0; i < L; i++) {
        Arr[i] *= K1;
        Arr[i] /= K2;
    }
}

// Funzione che moltiplica per una costante tutti gli elementi di un array a valori reali di lunghezza L
// K1 = numeratore
// K2 = denominatore
static void arr_scaling_re(float *Arr, int L, float K1, float K2) {
    for (int i = 0; i < L; i++) {
        Arr[i] *= K1;
        Arr[i] /= K2;
    }
}

static void mdf_adjust_prop(a3l_aec_state* st) {
    // anziché passare W, N e M, passo direttamente il puntatore alla struttura così può leggersi e modificarsi i campi senza restituire niente
    float tmp;
    for (int i = 0; i < st->M; i++) {
        st->prop[i] = 0.0f;                                                                
        tmp = 1.0f;
        for (int j = 0; j <= st->frame_size + 1; j++) {
            tmp += powf(cabsf(st->W[i][j]), 2.0f);
        }
        st->prop[i] = sqrtf(tmp);
    }
    st->prop_sum = 1.0f;
    for (int i = 0; i < st->M; i++) {
        st->max_sum[i] = 0.0f;
        st->max_sum[i] = (st->prop[i] > 1.0f) ? st->prop[i] : 1.0f;
        st->prop[i] += 0.1f*st->max_sum[i];
        st->prop_sum += st->prop[i];
        st->prop[i] = 0.99f*st->prop[i] / st->prop_sum;
    }
}

// Loop di gestione dell'elaborazione a blocchi
void a3l_aec_main_loop(a3l_aec_state *st, const int16_t *x, const int16_t *s, int16_t *e, fftwf_plan fft_plan, fftwf_plan ifft_plan) {
    // Parametri:
    // st: struttura contenente i parametri
    // x: segnale proveniente dal FE
    // s: segnale locale preso dal microfono
    // e: array di output

    float Pey_cur = 1.0f;                                                          //[a3l_mdf_1.m line 117]
    float Pyy_cur = 1.0f;                                                          //[a3l_mdf_1.m line 118]
    float Pey_tmp, Pyy_tmp;                                                        // nuovo nome delle variabili Pey e Pyy fuori dalla struttura)
    float Sxx, Sff, Dbf, See, Sey, Syy, Sdd;
    float ss, ss_1;
    float VAR1_UPDATE, VAR2_UPDATE, VAR_BACKTRACK, MIN_LEAK;
    float tmp_out, tmp32;
    float alpha, alpha_1;
    float RER;
    float r, er;
    float adapt_rate;

    for (int i = 0; i < st->frame_size; i++) {                                  
        st->out[i] = 0.0f;
    } 

    st->cancel_count++;                                                             //[a3l_mdf_1.m line 120]
                
    ss = 0.35f / st->M;                                                             //[a3l_mdf_1.m line 121]
    ss_1 = 1.0f - ss;                                                               //[a3l_mdf_1.m line 122]

    // memcpy(st->input, s_frame_f32, st->frame_size * sizeof(float));                 //[a3l_mdf_1.m line 123] override notch filtering

    for (int i = 0; i < st->frame_size; i++) {                                      //[a3l_mdf_1.m line 124-127]
        st->input[i] = (float)s[i];
        st->x[i] = st->x[i + st->frame_size];
        st->x[i + st->frame_size] = (float)x[i];
    }

    circshift(st);                                                                  //[a3l_mdf_1.m line 129]

    memcpy(st->ifft_buf, st->x, st->window_size * sizeof(float));
    fftwf_execute(fft_plan);                                                        //[a3l_mdf_1.m line 130]
    memcpy(st->X[0], st->fft_buf, st->window_size * sizeof(float _Complex));
    arr_scaling_cplx(st->X[0], st->window_size, 1.0f, (float)st->window_size);      // divide per 256
    conj_sym(st->X[0], st->window_size);                                            // rendo l'output coniugato simmetrico, poichè il piano fftw non lo fa

    Sxx = 0.0f;                                                                      //[a3l_mdf_1.m line 131]

    for (int i = 0; i < st->frame_size; i++) {                                      //[a3l_mdf_1.m line 132]
        Sxx = Sxx + powf(st->x[i + st->frame_size], 2.0f);
    }

    for (int i = 0; i < st->frame_size + 1; i++) {                                  //[a3l_mdf_1.m line 133]
        st->Xf[i] = powf(cabsf(st->X[0][i]), 2.0f);
    }

    Sff = 0.0f;                                                                      //[a3l_mdf_1.m line 131]
    
    for (int i = 0; i < st->window_size; i++) {                                     //[a3l_mdf_1.m line 135]
        st->Y[i] = 0.0f;
    }

    for (int j = 0; j < st->M; j++) {                                               //[a3l_mdf_1.m line 136-138]
        for (int i = 0; i < st->window_size; i++) {
            st->Y[i] += st->X[j][i] * st->foreground[j][i];
        }
    }

    memcpy(st->fft_buf, st->Y, st->window_size * sizeof(float _Complex));
    fftwf_execute(ifft_plan);                                                        //[a3l_mdf_1.m line 140]
    memcpy(st->e, st->ifft_buf, st->window_size * sizeof(float));

    for (int i = 0; i < st->frame_size; i++) {                                      //[a3l_mdf_1.m line 141]
        st->e[i] = st->input[i] - st->e[i + st->frame_size];
        Sff += powf(fabsf(st->e[i]), 2.0f);
    }

    if (st->adapted == 1){
        mdf_adjust_prop(st);
    }

    if (st->saturated == 0) {                                                       //[a3l_mdf_1.m line 148-155]
        float _Complex power_1_buf[256];
        for (int j = st->M - 1; j >= 0; j--) {
            for (int i = 1; i < st->frame_size; i++) {
                power_1_buf[i] = st->power_1[i];
                power_1_buf[st->window_size - i] = st->power_1[i];
            }
            power_1_buf[0] = st->power_1[0];
            power_1_buf[st->frame_size] = st->power_1[st->frame_size];

            for (int i = 0; i < st->window_size; i++) {
                st->PHI[i] = power_1_buf[i] * st->prop[j] * conjf(st->X[j + 1][i]) * st->E[i];
            }
            for (int i = 0; i < st->window_size; i++) {
                st->W[j][i] += st->PHI[i];
            }
        }
    } else {
        st->saturated -= 1;
    }

    for (int j = 0; j < st->M; j++) {                                               //[a3l_mdf_1.m line 157-163]
    #ifdef AUMDF                                                                    // compila se è definita la macro AUMDF
        if (j == 0 || (2 + st->cancel_count) % (st->M - 1) == j) {                  
            memcpy(st->fft_buf, st->W[j], st->window_size * sizeof(float _Complex));
            fftwf_execute(ifft_plan);
            memcpy(st->wtmp, st->ifft_buf, st->window_size * sizeof(float));
            arr_scaling_re(st->wtmp, st->window_size, 1.0f, (float)st->window_size);    // divido per 256
            for (int i = st->frame_size+1; i < st->window_size; i++) {
                st->wtmp[i] = 0.0f;
            }
            memcpy(st->ifft_buf, st->wtmp, st->window_size * sizeof(float));
            fftwf_execute(fft_plan);
            memcpy(st->W[j], st->fft_buf, st->window_size * sizeof(float _Complex));
            conj_sym(st->W[j], st->window_size);
        }
    #endif
    #ifdef MDF                                                                      // compila se è definita la macro MDF
        memcpy(st->fft_buf, st->W[j], st->window_size * sizeof(float _Complex));
        fftwf_execute(ifft_plan);
        memcpy(st->wtmp, st->ifft_buf, st->window_size * sizeof(float));
        arr_scaling_re(st->wtmp, st->window_size, 1.0f, (float)st->window_size);    // divido per 256
        for (int i = st->frame_size+1; i < st->window_size; i++) {
            st->wtmp[i] = 0.0f;
        }
        memcpy(st->ifft_buf, st->wtmp, st->window_size * sizeof(float));
        fftwf_execute(fft_plan);
        memcpy(st->W[j], st->fft_buf, st->window_size * sizeof(float _Complex));
        conj_sym(st->W[j], st->window_size);
    #endif
    }

    for (int i = 0; i < st->frame_size + 1; i++) {                                  //[a3l_mdf_1.m line 165-167]
        st->Yf[i] = 0.0f;
        st->Rf[i] = 0.0f;
        st->Xf[i] = 0.0f;
    }

    Dbf = 0.0f;                                                                     //[a3l_mdf_1.m line 168]
    
    for (int i = 0; i < st->window_size; i++) {                                     //[a3l_mdf_1.m line 169]
        st->Y[i] = 0.0f;
    }

    for (int j = 0; j < st->M; j++) {                                               //[a3l_mdf_1.m line 171-173]
        for (int i = 0; i < st->window_size; i++) {  
            st->Y[i] += st->X[j][i] * st->W[j][i];
        }
    }

    memcpy(st->fft_buf, st->Y, st->window_size * sizeof(float _Complex));           //[a3l_mdf_1.m line 175]
    fftwf_execute(ifft_plan);
    memcpy(st->y, st->ifft_buf, st->window_size * sizeof(float));

    See = 0.0f;                                                                      //[a3l_mdf_1.m line 176]

    for (int i = 0; i < st->frame_size; i++) {                                      
        st->e[i] = st->e[i + st->frame_size] - st->y[i + st->frame_size];            //[a3l_mdf_1.m line 177]
        Dbf += powf(fabsf(st->e[i]), 2.0f);                                          //[a3l_mdf_1.m line 178]
    }

    Dbf += 10.0f;

    for (int i = 0; i < st->frame_size; i++) {                                      
        st->e[i] = st->input[i] - st->y[i + st->frame_size];                        //[a3l_mdf_1.m line 179]
        See += powf(fabsf(st->e[i]), 2.0f);                                         //[a3l_mdf_1.m line 180]
    }

    VAR1_UPDATE = 0.5f;                                                              //[a3l_mdf_1.m line 181]
    VAR2_UPDATE = 0.25f;                                                             //[a3l_mdf_1.m line 182]
    VAR_BACKTRACK = 4.0f;                                                            //[a3l_mdf_1.m line 183]
    MIN_LEAK = 0.005f;                                                               //[a3l_mdf_1.m line 184]
    st->Davg1 = 0.6f*st->Davg1 + 0.4f*(Sff - See);                                   //[a3l_mdf_1.m line 185]
    st->Davg2 = 0.85f*st->Davg2 + 0.15f*(Sff - See);                                 //[a3l_mdf_1.m line 186]
    st->Dvar1 = 0.36f*st->Dvar1 + 0.16f*Sff*Dbf;                                     //[a3l_mdf_1.m line 187]
    st->Dvar2 = 0.7225f*st->Dvar2 + 0.0225f*Sff*Dbf;                                 //[a3l_mdf_1.m line 188]
    st->update_foreground = 0;                                                       //[a3l_mdf_1.m line 189]

    if ((Sff - See) * fabsf(Sff - See) > (Sff * Dbf)) {                             //[a3l_mdf_1.m line 190-197]        
        st->update_foreground = 1;
    } else if ((st->Davg1 * fabsf(st->Davg1)) > (VAR1_UPDATE * st->Dvar1)) {
        st->update_foreground = 1;
    } else if ((st->Davg2 * fabsf(st->Davg2)) > (VAR2_UPDATE * st->Dvar2)) {
        st->update_foreground = 1;
    }

    if (st->update_foreground == 1){                                                //[a3l_mdf_1.m line 199-206]
        st->Davg1 = 0.0f;
        st->Davg2 = 0.0f;
        st->Dvar1 = 0.0f;
        st->Dvar2 = 0.0f;
        for (int i = 0; i < st->M; i++) {
            memcpy(st->foreground[i], st->W[i], st->window_size * sizeof(float _Complex));
        }
        for (int i = 0; i < st->frame_size; i++) {                                      
            st->e[st->frame_size + i] = (st->window[st->frame_size + i] * st->e[st->frame_size + i]) + (st->window[i] * st->y[st->frame_size + i]);                                                            
        }
    } else {
        st->reset_background = 0;                                                   //[a3l_mdf_1.m line 207]

        if (-(Sff-See)*fabsf(Sff-See) > VAR_BACKTRACK*(Sff*Dbf)){                   //[a3l_mdf_1.m line 209-211]
            st->reset_background = 1;
        }

        if ((-st->Davg1*fabsf(st->Davg1)) > VAR_BACKTRACK*st->Dvar1){               //[a3l_mdf_1.m line 213-215]
            st->reset_background = 1;
        }

        if ((-st->Davg2*fabsf(st->Davg2)) > VAR_BACKTRACK*st->Dvar2){               //[a3l_mdf_1.m line 217-219]
            st->reset_background = 1;
        }

        if (st->reset_background == 1){
            for (int i = 0; i < st->M; i++) {                                       //[a3l_mdf_1.m line 222]
                memcpy(st->W[i], st->foreground[i], st->window_size * sizeof(float _Complex));
            }
            for (int i = 0; i < st->frame_size; i++) {                              
                st->y[st->frame_size + i] = st->e[st->frame_size + i];              //[a3l_mdf_1.m line 223]
                st->e[i] = st->input[i] - st->y[st->frame_size + i];                //[a3l_mdf_1.m line 224]
            }
            See = Sff;                                                              //[a3l_mdf_1.m line 225-229]
            st->Davg1 = 0.0f;
            st->Davg2 = 0.0f;
            st->Dvar1 = 0.0f;
            st->Dvar2 = 0.0f;
        }
    }

    Sey = 0.0f;                                                                     //[a3l_mdf_1.m line 233-235]
    Syy = 0.0f;
    Sdd = 0.0f;

    for (int i = 0; i < st->frame_size; i++) {                                      //[a3l_mdf_1.m line 237-246]
        tmp_out = st->input[i] - st->e[st->frame_size + i];
        if (st->input[i] <= -32000.0f || st->input[i] >= 32000.0f){
            if (st->saturated == 0){
                st->saturated = 1;
            }
        }
        st->out[i] = tmp_out;
    }

    for (int i = 0; i < st->frame_size; i++) {
        st->e[st->frame_size + i] = st->e[i];                                       //[a3l_mdf_1.m line 248]
        st->e[i] = 0.0f;                                                             //[a3l_mdf_1.m line 249]
        Sey += (st->e[st->frame_size + i] * st->y[st->frame_size + i]);             //[a3l_mdf_1.m line 250]
        Syy += powf(st->y[st->frame_size + i], 2.0f);                                   //[a3l_mdf_1.m line 251]
        Sdd += powf(st->input[i], 2.0f);                                                //[a3l_mdf_1.m line 252]
    }

    memcpy(st->ifft_buf, st->e, st->window_size * sizeof(float));
    fftwf_execute(fft_plan);                                                         //[a3l_mdf_1.m line 253]
    memcpy(st->E, st->fft_buf, st->window_size * sizeof(float _Complex));
    arr_scaling_cplx(st->E, st->window_size, 1.0f, (float)st->window_size);         // divide per 256
    conj_sym(st->E, st->window_size);

    for (int i = 0; i < st->frame_size; i++) {
        st->y[i] = 0.0f;                                                             //[a3l_mdf_1.m line 254]
    }

    memcpy(st->ifft_buf, st->y, st->window_size * sizeof(float));
    fftwf_execute(fft_plan);                                                         //[a3l_mdf_1.m line 255]
    memcpy(st->Y, st->fft_buf, st->window_size * sizeof(float _Complex));
    arr_scaling_cplx(st->Y, st->window_size, 1.0f, (float)st->window_size);
    conj_sym(st->Y, st->window_size);

    for (int i = 0; i < st->frame_size + 1; i++) {                                  //[a3l_mdf_1.m line 256-256]
        st->Rf[i] += powf(cabsf(st->E[i]), 2.0f);
        st->Yf[i] += powf(cabsf(st->Y[i]), 2.0f);
    }

    if (!(Syy >= 0 && Sxx >= 0 && See >= 0)){                                       //[a3l_mdf_1.m line 259-266]
        st->screwed_up += 50;
        for (int i = 0; i < st->frame_size; i++) {                                  
            st->out[i] = 0.0;
        } 
    } else if (Sff > Sdd + (float)(st->window_size*10000.0f)) {
        st->screwed_up += 1;
    } else {
        st->screwed_up = 0;
    }

    if (st->screwed_up >= 50){                                                      //[a3l_mdf_1.m line 268-271]
        printf(">> AEC full reset!\n");
        a3l_aec_state_reset(st);
    }

    See = fmax(See, (float)(st->window_size*100.0f));
    
    for (int i = 0; i < st->frame_size; i++) {
        Sxx += powf(st->x[st->frame_size + i], 2.0f);                                   //[a3l_mdf_1.m line 274]
    }

    for (int i = 0; i < st->frame_size + 1; i++) {
        st->Xf[i] = powf(cabsf(st->X[0][i]), 2.0f);                                      //[a3l_mdf_1.m line 275]
        st->power[i] = ss_1*st->power[i] + 1.0f + ss*st->Xf[i];                      //[a3l_mdf_1.m line 276]
        st->Eh_cur[i] = st->Rf[i] - st->Eh[i];                                      //[a3l_mdf_1.m line 277]
        st->Yh_cur[i] = st->Yf[i] - st->Yh[i];                                      //[a3l_mdf_1.m line 278]
        Pey_cur += st->Eh_cur[i]*st->Yh_cur[i];                                     //[a3l_mdf_1.m line 279]
        Pyy_cur += st->Yh_cur[i]*st->Yh_cur[i];                                     //[a3l_mdf_1.m line 280]
        st->Eh[i] = (1.0f - st->spec_average)*st->Eh[i] + st->spec_average*st->Rf[i];//[a3l_mdf_1.m line 281]
        st->Yh[i] = (1.0f - st->spec_average)*st->Yh[i] + st->spec_average*st->Yf[i];//[a3l_mdf_1.m line 282]
    }

    Pyy_tmp = sqrtf(Pyy_cur);                                                        //[a3l_mdf_1.m line 283]
    Pey_tmp = Pey_cur / Pyy_tmp;                                                    //[a3l_mdf_1.m line 284]
    tmp32 = st->beta0 * Syy;                                                        //[a3l_mdf_1.m line 285]

    if (tmp32 > st->beta_max*See){                                                  //[a3l_mdf_1.m line 287-289]
        tmp32 = st->beta_max*See;
    }

    alpha = tmp32 / See;                                                            //[a3l_mdf_1.m line 291]
    alpha_1 = 1.0f - alpha;                                                          //[a3l_mdf_1.m line 292]
    st->Pey = alpha_1*st->Pey + alpha*Pey_tmp;                                      //[a3l_mdf_1.m line 293]
    st->Pyy = alpha_1*st->Pyy + alpha*Pyy_tmp;                                      //[a3l_mdf_1.m line 294]

    if (st->Pyy < 1.0f){                                                             //[a3l_mdf_1.m line 296-298]
        st->Pyy = 1.0f;
    }

    if (st->Pey < MIN_LEAK*st->Pyy){                                                //[a3l_mdf_1.m line 300-302]
        st->Pey = MIN_LEAK*st->Pyy;
    }

    if (st->Pey > st->Pyy){                                                         //[a3l_mdf_1.m line 304-306]
        st->Pey = st->Pyy;
    }

    st->leak_estimate = st->Pey / st->Pyy;                                          //[a3l_mdf_1.m line 308]

    if (st->leak_estimate > 16383.0f){                                                 //[a3l_mdf_1.m line 310-312]
        st->leak_estimate = 32767.0f;
    }

    RER = (0.0001f*Sxx + 3.0f*st->leak_estimate*Syy) / See;                           //[a3l_mdf_1.m line 314]

    if (RER < (Sey*Sey/(1+See*Syy))){                                               //[a3l_mdf_1.m line 316-318]
        RER = (Sey*Sey/(1+See*Syy));
    }

    if (RER > 0.5f){                                                                 //[a3l_mdf_1.m line 320-322]
        RER = 0.5f;
    }

    if (st->adapted==0 && st->sum_adapt > st->M && st->leak_estimate * Syy > 0.03f * Syy) {//[a3l_mdf_1.m line 324-326]
        st->adapted = 1;
    }

    if (st->adapted == 1){                                                          //[a3l_mdf_1.m line 328-349]
        for (int i = 0; i < st->frame_size + 1; i++) {
            r = st->leak_estimate * st->Yf[i];
            er = st->Rf[i] + 1.0f;
            if (r > 0.5f*er){
                r = 0.5f*er;
            }
            r = 0.7f*r + 0.3f*RER*er;
            st->power_1[i] = r / (er*st->power[i]+10.0f);
        }
    } else {                    
        adapt_rate = 0.0f;
        if (Sxx > st->window_size*1000) {
            tmp32 = 0.25f * Sxx;
            if (tmp32 > 0.25f*See) {
                tmp32 = 0.25f*See;
            }
            adapt_rate = tmp32 / See;
        }
        for (int i = 0; i < st->frame_size + 1; i++) {                              //[a3l_mdf_1.m line 347]
            st->power_1[i] = adapt_rate / (st->power[i] + 10.0f);
        }
        st->sum_adapt += adapt_rate;
    }

    for (int i = 0; i < st->frame_size; i++) {                                      //[a3l_mdf_1.m line 351]
        st->last_y[i] = st->last_y[st->frame_size + i];
    }

    if (st->adapted == 1) {                                                         //[a3l_mdf_1.m line 353-355]
        for (int i = 0; i < st->frame_size; i++) {                                  
            st->last_y[st->frame_size + i] = st->input[i] - st->out[i];
        }
    }
    
    // clip and cast to short
    for (int i = 0; i < st->frame_size; i++) {
		if (st->out[i] > 32767.0f) {
			e[i] = 32767;
		} else if (st->out[i] < -32768.0f) {
			e[i] = -32768;
		} else {
			e[i] = (int16_t)(st->out[i]);
		}
    }
}

