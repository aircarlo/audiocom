// a3l_enhancement.h
#ifndef ENH_H
#define ENH_H

#include <stdint.h>

// Definizione della struttura dati per le funzioni di noise cancellation
typedef struct {                        
	int frame_size;                     // 128
	int window_size;                    // 256
	float vad_decision;
	float *noise_PSD_buf;
	float *Y_psd_buf;
	float *yf;
    float *yf_buf;
	float *out_buf;
	float *xf_hat;
    float *hann_win;
	float _Complex *Yc;
	float *Y_mag;
	float *Y_psd;
	float *gamma_k;
	float *csi_k;
	float *log_sigma_k;
	float *Y_pha;
	float _Complex *X_hat;
    float *A;
    float *vk;
    float *ei_vk;
    float *hw;
    float _Complex *fft_buf;
    float *ifft_buf;
} a3l_se_state;

// Funzione di allocazione dinamica e inizializzazione
void a3l_se_state_init(a3l_se_state *st, int frame_size);

// funzione di deallocazione
void a3l_se_state_free(a3l_se_state *st);

// Funzione di elaborazione principale
void a3l_enh_main_loop(a3l_se_state* st, const int16_t *in, int16_t *out, fftwf_plan fft_plan, fftwf_plan ifft_plan);

#endif // ENH_H