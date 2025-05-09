// a3l_aec.h
#ifndef AEC_H
#define AEC_H

#include <stdint.h>

 // Struttura per i parametri dell'algoritmo di cancellazione
 typedef struct {
    int frame_size;
    int window_size;
    int M;
    int update_foreground;              // usato come flag 0/1. In Matlab non faceva parte di st
    int reset_background;               // usato come flag 0/1. In Matlab non faceva parte di st
    int cancel_count;
    int saturated;
    int screwed_up;
    int adapted;
    float sum_adapt;
    float spec_average;
    float beta0;
    float beta_max;
    float leak_estimate;
    float prop_sum;                    // variabile interna alla funzione mdf_adjust_prop
    float *e;
    float *x;
    float *input;
    float *out;
    float *y;
    float *last_y;
    float *Yf;
    float *Rf;
    float *Xf;
    float *Yh;
    float *Eh;
    float *Eh_cur;
    float *Yh_cur;
    float *ifft_buf;                   // buffer di appoggio dopo la ifft
    float _Complex *fft_buf;           // buffer di appoggio dopo la fft
    float _Complex **X;
    float _Complex *Y;
    float _Complex *E;
    float _Complex **W;
    float _Complex **foreground;
    float _Complex *PHI;
    float *power;
    float *power_1;
    float *window;
    float *prop;
    float *max_sum;
    float *wtmp;
    float Pey;
    float Pyy;
    float Davg1;
    float Davg2;
    float Dvar1;
    float Dvar2;
} a3l_aec_state;

// Funzione di allocazione dinamica e inizializzazione
void a3l_aec_state_init(a3l_aec_state *st, int frame_size, int filter_length, int sample_rate);

// Funzione di reset della struttura (in caso di fault)
void a3l_aec_state_reset(a3l_aec_state *st);

// funzione di deallocazione
void a3l_aec_state_free(a3l_aec_state *st);

// Funzione di elaborazione principale
void a3l_aec_main_loop(a3l_aec_state *st, const int16_t *x, const int16_t *s, int16_t *e, fftwf_plan fft_plan, fftwf_plan ifft_plan);

#endif // AEC_H